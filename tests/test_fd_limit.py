#!/usr/bin/env python3
"""Regression test for the dante_ltr "Too many open files" chunk-split bug.

The chunk split (main(), the `[open(f) for f in temp_files]` loops) opens one
file handle per chunk *simultaneously*. The number of chunks is
`int(total_size / max_chunk_size)`, so on large genomes it exceeds the open-file
limit and the run aborts with `OSError: [Errno 24] Too many open files`.

`_limit_chunks_to_fd_budget()` must keep the chunk count within the (raised) fd
limit so the split always fits. This test forces a low limit and checks that:
  1. the guard caps the chunk count below the limit,
  2. that capped number of handles can actually be opened at once, and
  3. the original un-capped count would indeed have exhausted the fds.
It also checks the guard is a no-op when the limit is ample.
"""
import importlib.machinery
import os
import resource
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)  # so the script's `from version import __version__` resolves


def load_dante_ltr():
    loader = importlib.machinery.SourceFileLoader(
        "dante_ltr_mod", os.path.join(ROOT, "dante_ltr")
    )
    return loader.load_module()


def open_all(n):
    """Open n fresh temp files at once (the original crashing pattern)."""
    handles, paths = [], []
    try:
        for _ in range(n):
            fd, path = tempfile.mkstemp()
            os.close(fd)
            paths.append(path)
            handles.append(open(path, "w"))
        return True
    finally:
        for h in handles:
            h.close()
        for p in paths:
            os.remove(p)


def run_constrained(limit):
    # Lower BOTH soft and hard so the guard cannot raise its way out and must cap.
    resource.setrlimit(resource.RLIMIT_NOFILE, (limit, limit))
    mod = load_dante_ltr()
    requested = limit * 4  # far more chunks than the fd budget
    allowed = mod._limit_chunks_to_fd_budget(requested)
    assert allowed < limit, f"guard did not cap: {allowed} vs limit {limit}"
    assert open_all(allowed), "opening the capped number of handles failed"
    crashed = False
    try:
        open_all(requested)
    except OSError:
        crashed = True
    assert crashed, "expected the un-capped count to exhaust fds under the low limit"
    print(f"  constrained(limit={limit}): requested {requested} -> capped {allowed}, opens OK")


def run_ample():
    mod = load_dante_ltr()
    n = 50
    allowed = mod._limit_chunks_to_fd_budget(n)
    assert allowed == n, f"ample case changed the count: {allowed} != {n}"
    print(f"  ample: {n} -> {allowed} (unchanged)")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "--child":
        run_constrained(int(sys.argv[2]))
        sys.exit(0)
    run_ample()
    # Run the constrained case in a child so lowering the hard limit (irreversible
    # within a process) does not affect this test runner.
    import subprocess

    subprocess.run(
        [sys.executable, os.path.abspath(__file__), "--child", "128"], check=True
    )
    print("test_fd_limit: PASSED")
