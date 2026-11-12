# Launch the web app

PathwayAgent runs **on your own computer**. This page checks whether it is running and opens it.
Your genes are not uploaded to this website.

<div class="pa-launch" markdown>
<p id="pa-status" class="pa-status pa-checking">Checking for PathwayAgent on this computer …</p>

[Open PathwayAgent](http://127.0.0.1:8765){ #pa-open .md-button .md-button--primary target=_blank }
<button id="pa-recheck" class="md-button" type="button">Check again</button>
</div>

In the app, paste genes or upload a file (gene list, ranked list, DE table, or one module per
column), choose the organism and preset, and run. See
[Inputs and command line](usage.md#input-formats) for the formats and
[Web app and REST API](webapp.md) for what the page offers.

## Start PathwayAgent

If the status above says it is not running, start it in one of two ways. Both need the repository
on your computer first ([download or clone it](https://github.com/helenhuangmath/PathwayAgent)).

=== "Double-click (macOS)"

    Double-click **`PathwayAgent.command`** in the repository folder.

    The first start creates a Python environment and builds the gene-set library (about five
    minutes, needs internet and Python 3.11 or newer). Later starts take a few seconds. A Terminal
    window stays open while the app runs; close it to stop.

    If macOS refuses to open the file because it was downloaded from the internet, right-click it,
    choose **Open**, and confirm once.

=== "Terminal (macOS, Linux)"

    After [installing](INSTALL.md):

    ```bash
    .venv/bin/python -m pathwayagent serve
    ```

    or, without installing first, let the launcher set everything up:

    ```bash
    bash PathwayAgent.command
    ```

The app opens at <http://127.0.0.1:8765>. Then come back here and press **Check again**, or use
the address directly.

## Good to know

- **Privacy.** The app listens on `127.0.0.1` only, so nobody else can reach it. Tools that use a
  web service (Enrichr, g:Profiler, Reactome, STRING) still send your genes to that service unless
  you tick **local only**.
- **The status check.** Some browsers do not let a website look for programs on your computer, or
  ask for permission first ("access devices on your local network"). If the status stays at "not
  detected" while the app is running, **Open PathwayAgent** still works.
- **Another port.** If you started the app with `--port`, open `http://127.0.0.1:<port>` yourself.
- **Peaks.** Genomic regions (BED) are not accepted yet; convert peaks to genes first.

<script>
(function () {
  var url = "http://127.0.0.1:8765";
  var status = document.getElementById("pa-status");
  var recheck = document.getElementById("pa-recheck");
  if (!status || !recheck) return;

  function show(kind, text) {
    status.className = "pa-status pa-" + kind;
    status.textContent = text;
  }

  function check() {
    show("checking", "Checking for PathwayAgent on this computer …");
    var ctl = new AbortController();
    var timer = setTimeout(function () { ctl.abort(); }, 3000);
    fetch(url + "/api/v1/health", { signal: ctl.signal, cache: "no-store" })
      .then(function (r) { return r.json(); })
      .then(function (j) {
        if (!j || !j.ok) throw new Error("unexpected reply");
        show("up", "PathwayAgent " + (j.version || "") + " is running on this computer. Press Open PathwayAgent.");
      })
      .catch(function () {
        show("down", "PathwayAgent was not detected on this computer. Start it as described below, then check again.");
      })
      .then(function () { clearTimeout(timer); });
  }

  recheck.addEventListener("click", check);
  check();
})();
</script>
