#!/bin/bash
# PathwayAgent launcher. Double-click this file on macOS (or run `bash PathwayAgent.command` on
# Linux) to start the local web app and open it in the browser. On first use it creates a Python
# environment in .venv/ and builds the gene-set library; later starts take a few seconds.
#
# Environment variables: PATHWAYAGENT_PORT (default 8765), PATHWAYAGENT_PYTHON (Python 3.11+ to
# create the environment with), and the PATHWAYAGENT_* paths described in docs/INSTALL.md.

cd "$(dirname "$0")" || exit 1
PORT="${PATHWAYAGENT_PORT:-8765}"
URL="http://127.0.0.1:${PORT}"

pause_and_exit() {
    echo
    read -r -p "Press Return to close this window. " _
    exit "${1:-1}"
}

open_url() {
    if command -v open >/dev/null 2>&1; then open "$1"
    elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$1" >/dev/null 2>&1
    else echo "Open $1 in your browser."
    fi
}

echo "PathwayAgent"
echo "============"

# Already running: just open it.
if curl -fs --max-time 2 "${URL}/api/v1/health" 2>/dev/null | grep -q '"ok"'; then
    echo "PathwayAgent is already running at ${URL}."
    open_url "${URL}"
    exit 0
fi

# First use: create the environment and install the package.
if [ ! -x .venv/bin/python ]; then
    PY=""
    for cand in "${PATHWAYAGENT_PYTHON}" python3.12 python3.11 python3.13 python3; do
        [ -n "${cand}" ] || continue
        if command -v "${cand}" >/dev/null 2>&1 &&
           "${cand}" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
            PY="${cand}"
            break
        fi
    done
    if [ -z "${PY}" ]; then
        echo "Python 3.11 or newer was not found. Install it from https://www.python.org/downloads/"
        echo "and start this launcher again."
        pause_and_exit 1
    fi
    echo "First start: installing PathwayAgent into .venv/ with $(${PY} --version) (a few minutes) ..."
    "${PY}" -m venv .venv || pause_and_exit 1
    if ! .venv/bin/python -m pip install --quiet --upgrade pip ||
       ! .venv/bin/python -m pip install --quiet -e . scikit-learn; then
        echo "Installation failed. See the messages above."
        rm -rf .venv
        pause_and_exit 1
    fi
fi

# First use: build the gene-set library.
if ! .venv/bin/python -m pathwayagent db list 2>/dev/null | grep -q ": built "; then
    echo "Building the gene-set library (one time, about 5 minutes, needs internet) ..."
    .venv/bin/python -m pathwayagent db build || { echo "The library build failed."; pause_and_exit 1; }
fi

echo "Starting PathwayAgent at ${URL}"
echo "Keep this window open while you work. Close it, or press Ctrl+C, to stop."
echo
.venv/bin/python -m pathwayagent serve --port "${PORT}"
status=$?
[ "${status}" -eq 0 ] || [ "${status}" -eq 130 ] || pause_and_exit "${status}"
