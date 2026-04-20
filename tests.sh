#!/bin/bash
# tests.sh — dispatcher. Run a specific test level:
#   ./tests.sh smoke         # < 30 s, runs on every PR
#   ./tests.sh short         # ~1 min, runs on every PR
#   ./tests.sh fallback      # ~1-2 min, runs on every PR
#   ./tests.sh long          # ~10-30 min, runs on release tags
#   ./tests.sh all           # smoke + short + fallback + long
#
# Backwards compat: if the first argument is a number, treat it as CPU
# count and run the long test (old behaviour of ./tests.sh 4).
#
# NCPU can also be set via the NCPU env var.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

LEVEL="${1:-smoke}"
export NCPU="${NCPU:-${2:-2}}"

# Numeric first arg -> legacy mode (long test)
if [[ "$LEVEL" =~ ^[0-9]+$ ]]; then
  export NCPU="$LEVEL"
  LEVEL="long"
fi

# Some environments need an explicit conda activate; but in CI we run
# inside an already-activated environment. Activate only if conda is on PATH
# and we are not already in the env.
if command -v conda >/dev/null 2>&1 && [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
  eval "$(conda shell.bash hook)"
  if conda env list | grep -q '^dante_ltr '; then
    conda activate dante_ltr
  fi
fi

case "$LEVEL" in
  smoke)    bash "$ROOT/tests/smoke.sh" ;;
  short)    bash "$ROOT/tests/short.sh" ;;
  fallback) bash "$ROOT/tests/fallback.sh" ;;
  long)     bash "$ROOT/tests/long.sh"  ;;
  all)
    bash "$ROOT/tests/smoke.sh"
    bash "$ROOT/tests/short.sh"
    bash "$ROOT/tests/fallback.sh"
    bash "$ROOT/tests/long.sh"
    ;;
  *)
    echo "usage: $0 {smoke|short|fallback|long|all|<NCPU>}" >&2
    exit 2
    ;;
esac
