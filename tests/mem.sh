#!/bin/bash
# tests/mem.sh — regression test for the chunk-pool memory-gate budget (issue #13).
# The gate sized the pool from /proc/meminfo MemAvailable, which the kernel does not
# namespace: inside a container under a batch scheduler it reports the whole node, so
# a memory-limited job got a pool many times too large and was OOM-killed. Verify the
# budget resolution chain (--max_memory / scheduler env / cgroup / MemAvailable) and
# the host-wide-budget warning. Fast (<1 s), no genome needed.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
echo "=== memory budget resolution (chunk pool OOM regression) ==="
python3 tests/test_mem_budget.py
echo
echo "mem PASSED"
