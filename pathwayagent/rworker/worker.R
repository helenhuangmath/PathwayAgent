# PathwayAgent persistent R worker (design §7.3).
# Protocol: one JSON request per stdin line  {"id":..., "task":..., "args":{...}}
#           one JSON reply per request on stdout, prefixed with @@JSON@@
# Packages are loaded once, gene-set GMTs are cached, so each call costs only the test itself.

suppressPackageStartupMessages({
  library(jsonlite)
  library(clusterProfiler)
  library(fgsea)
})
options(warn = -1)

emit <- function(obj) {
  cat("@@JSON@@", toJSON(obj, auto_unbox = TRUE, dataframe = "rows", null = "null",
                         na = "null", digits = NA), "\n", sep = "")
  flush(stdout())
}

.gmt <- new.env()
read_gmt <- function(path) {
  if (!exists(path, envir = .gmt, inherits = FALSE)) {
    parts <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
    ids <- vapply(parts, `[`, "", 1)
    genes <- lapply(parts, function(x) x[-(1:2)])
    t2g <- data.frame(term = rep(ids, lengths(genes)), gene = unlist(genes), stringsAsFactors = FALSE)
    t2n <- data.frame(term = ids, name = vapply(parts, `[`, "", 2), stringsAsFactors = FALSE)
    assign(path, list(t2g = t2g, t2n = t2n, sets = setNames(genes, ids)), envir = .gmt)
  }
  get(path, envir = .gmt, inherits = FALSE)
}

orgdb <- function(org) {
  pkg <- if (org == "mmusculus") "org.Mm.eg.db" else "org.Hs.eg.db"
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  get(pkg)
}

as_rows <- function(res, db) {
  df <- as.data.frame(res)
  if (is.null(df) || nrow(df) == 0) return(list())
  df$database <- db
  df
}

task_enricher <- function(a) {
  out <- list()
  uni <- if (length(a$universe)) unlist(a$universe) else NULL
  for (db in names(a$gmts)) {
    g <- read_gmt(a$gmts[[db]])
    res <- enricher(unlist(a$genes), universe = uni, TERM2GENE = g$t2g, TERM2NAME = g$t2n,
                    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                    minGSSize = a$min_size, maxGSSize = a$max_size)
    if (!is.null(res)) out[[db]] <- as_rows(res, db)
  }
  out
}

# Native GO annotation from org.*.eg.db (GOALL = propagated, all evidence codes, as enrichGO uses).
# enrichGO() rebuilds this table on every call (~65 s); we build it once per worker and reuse it.
.godata <- new.env()
go_tables <- function(org, ont, cache_dir = NULL) {
  key <- paste(org, ont)
  rds <- if (length(cache_dir)) file.path(cache_dir, paste0("orgdb_go_", ont, ".rds")) else ""
  if (!exists(key, envir = .godata, inherits = FALSE) && nzchar(rds) && file.exists(rds)) {
    assign(key, readRDS(rds), envir = .godata)
  }
  if (!exists(key, envir = .godata, inherits = FALSE)) {
    suppressPackageStartupMessages(library(GO.db))
    db <- orgdb(org)
    m <- suppressMessages(AnnotationDbi::select(db, keys = AnnotationDbi::keys(db, "GOALL"), keytype = "GOALL",
                                                columns = c("SYMBOL", "ONTOLOGYALL")))
    m <- unique(m[m$ONTOLOGYALL == ont & !is.na(m$SYMBOL), c("GOALL", "SYMBOL")])
    names(m) <- c("term", "gene")
    ids <- unique(m$term)
    t2n <- data.frame(term = ids, name = unname(AnnotationDbi::Term(ids)), stringsAsFactors = FALSE)
    assign(key, list(t2g = m, t2n = t2n), envir = .godata)
    if (nzchar(rds)) { dir.create(dirname(rds), showWarnings = FALSE, recursive = TRUE); saveRDS(get(key, envir = .godata), rds) }
  }
  get(key, envir = .godata, inherits = FALSE)
}

task_enrichgo <- function(a) {
  out <- list()
  uni <- if (length(a$universe)) unlist(a$universe) else NULL
  for (ont in unlist(a$onts)) {
    g <- go_tables(a$organism, ont, a$cache_dir)
    res <- enricher(unlist(a$genes), universe = uni, TERM2GENE = g$t2g, TERM2NAME = g$t2n,
                    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                    minGSSize = a$min_size, maxGSSize = a$max_size)
    if (!is.null(res)) out[[paste0("GO_", ont)]] <- as_rows(res, paste0("GO_", ont))
  }
  out
}

task_reactomepa <- function(a) {
  suppressPackageStartupMessages(library(ReactomePA))
  db <- orgdb(a$organism)
  ids <- suppressMessages(bitr(unlist(a$genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = db))
  uni <- NULL
  if (length(a$universe)) {
    uni <- suppressMessages(bitr(unlist(a$universe), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = db))$ENTREZID
  }
  res <- enrichPathway(ids$ENTREZID, organism = if (a$organism == "mmusculus") "mouse" else "human",
                       universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                       minGSSize = a$min_size, maxGSSize = a$max_size, readable = TRUE)
  list(REACTOME = if (is.null(res)) list() else as_rows(res, "REACTOME"))
}

task_fgsea <- function(a) {
  stats <- setNames(as.numeric(unlist(a$scores)), unlist(a$genes))
  stats <- stats[!duplicated(names(stats))]
  out <- list()
  for (db in names(a$gmts)) {
    g <- read_gmt(a$gmts[[db]])
    set.seed(1)
    res <- fgseaMultilevel(g$sets, stats, minSize = a$min_size, maxSize = a$max_size, eps = 0)
    if (nrow(res)) {
      res$leadingEdge <- vapply(res$leadingEdge, paste, "", collapse = ";")
      res$database <- db
      out[[db]] <- as.data.frame(res)
    }
  }
  out
}

emit(list(ready = TRUE, versions = list(
  R = paste(R.version$major, R.version$minor, sep = "."),
  clusterProfiler = as.character(packageVersion("clusterProfiler")),
  fgsea = as.character(packageVersion("fgsea")))))

con <- file("stdin", open = "r")
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if (!nzchar(line)) next
  req <- tryCatch(fromJSON(line, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(req)) { emit(list(ok = FALSE, error = "bad json")); next }
  t0 <- Sys.time()
  res <- tryCatch({
    f <- switch(req$task, ping = function(a) list(pong = TRUE), enricher = task_enricher,
                enrichgo = task_enrichgo, reactomepa = task_reactomepa, fgsea = task_fgsea,
                stop(paste("unknown task", req$task)))
    list(ok = TRUE, result = f(req$args))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  res$id <- req$id
  res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  emit(res)
}
