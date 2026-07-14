#!/bin/bash
# tests/fd.sh — regression test for the "Too many open files" chunk-split bug.
# The chunk split opened one handle per chunk at once; on large genomes that
# count (total_size / max_chunk_size) exceeds ulimit -n and the run aborted with
# OSError: [Errno 24]. Verify the fd-budget guard keeps the split within the
# limit. Fast (<1 s), no genome needed.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
echo "=== fd budget guard (Too many open files regression) ==="
python3 tests/test_fd_limit.py
echo
echo "fd PASSED"
