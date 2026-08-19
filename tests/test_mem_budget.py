#!/usr/bin/env python3
"""Regression test for the dante_ltr chunk-pool memory-gate budget (issue #13).

The pool gate reads a memory budget and divides it by the measured per-chunk peak
RSS. It used to read `/proc/meminfo` `MemAvailable` only -- which the kernel does
not namespace, so inside a container under a batch scheduler it reports the whole
node rather than the job's allocation. `pool_size` is linear in that number, so a
128 GB job on a 1.6 TB node got a pool sized ~12x too large and was OOM-killed
hours in.

`_memory_budget_kb()` must resolve the budget from `--max_memory`, then the
scheduler environment, then the cgroup limit, and only then fall back to
`MemAvailable`. This test drives the helper against a synthetic /sys/fs/cgroup
tree and injected environments, checking each rung in isolation, the precedence
between them, and that `_warn_if_host_wide_budget()` fires only when the budget
really is host-wide *and* a scheduler is in charge.

The real inputs that exposed this (a 94 Gbp genome, 1889 chunks) cannot run on an
ordinary runner, which is why the helper takes injectable roots at all.
"""
import importlib.machinery
import io
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)  # so the script's `from version import __version__` resolves

KB_PER_GB = 1024 * 1024


def load_dante_ltr():
    loader = importlib.machinery.SourceFileLoader(
        "dante_ltr_mod", os.path.join(ROOT, "dante_ltr")
    )
    return loader.load_module()


mod = load_dante_ltr()


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def meminfo_file(tmp, kb):
    path = os.path.join(tmp, "meminfo")
    write(path, "MemTotal:  99999999 kB\nMemAvailable: {} kB\nCached: 1 kB\n".format(kb))
    return path


def budget(tmp, explicit=None, environ=None, cgroup_tree=None, proc_cgroup=None,
           meminfo_kb=None):
    """Call the helper with every filesystem root pointed at the synthetic tree.

    Roots that the case does not exercise are pointed at a path that cannot exist,
    so a rung that should miss cannot accidentally read the real host.
    """
    missing = os.path.join(tmp, "does-not-exist")
    return mod._memory_budget_kb(
        explicit_gb=explicit,
        environ={} if environ is None else environ,
        sysfs_root=cgroup_tree or missing,
        proc_cgroup=proc_cgroup or missing,
        meminfo=meminfo_file(tmp, meminfo_kb) if meminfo_kb else missing,
    )


def approx(a, b, tol=1):
    return abs(a - b) <= tol


def check(name, got, want_kb, want_source):
    kb, source = got
    assert source == want_source, "{}: source {!r} != {!r}".format(name, source, want_source)
    if want_kb is None:
        assert kb is None, "{}: expected no budget, got {}".format(name, kb)
    else:
        assert kb is not None and approx(kb, want_kb), \
            "{}: budget {} != {}".format(name, kb, want_kb)
    print("  {}: {} kB from {}".format(name, kb, source))


def test_explicit_flag_wins(tmp):
    # Every other rung is also satisfiable here; the flag must still win.
    env = {"AGENT_MEMORY": "64", "PBS_RESC_MEM": str(200 * 1024**3),
           "SLURM_MEM_PER_NODE": "300000"}
    check("explicit --max_memory", budget(tmp, explicit=128, environ=env,
                                          meminfo_kb=1600 * KB_PER_GB),
          128 * KB_PER_GB, "--max_memory")
    # Fractional GB, and the flag beating a host-wide MemAvailable.
    check("explicit fractional", budget(tmp, explicit=1.5), int(1.5 * KB_PER_GB),
          "--max_memory")


def test_agent_memory(tmp):
    check("AGENT_MEMORY", budget(tmp, environ={"AGENT_MEMORY": "16"}),
          16 * KB_PER_GB, "AGENT_MEMORY")


def test_scheduler_vars(tmp):
    # Each variable in its own unit, alone -- a 1024x unit slip here re-creates the bug.
    check("PBS_RESC_MEM (bytes)",
          budget(tmp, environ={"PBS_RESC_MEM": str(128 * 1024**3)}),
          128 * KB_PER_GB, "PBS_RESC_MEM")
    check("SLURM_MEM_PER_NODE (MB)",
          budget(tmp, environ={"SLURM_MEM_PER_NODE": "65536"}),
          64 * KB_PER_GB, "SLURM_MEM_PER_NODE")
    check("LSB_MAX_MEM_RUSAGE (kB)",
          budget(tmp, environ={"LSB_MAX_MEM_RUSAGE": str(32 * KB_PER_GB)}),
          32 * KB_PER_GB, "LSB_MAX_MEM_RUSAGE")
    check("SLURM_MEM_PER_CPU x CPUS_ON_NODE",
          budget(tmp, environ={"SLURM_MEM_PER_CPU": "4096",
                               "SLURM_CPUS_ON_NODE": "16"}),
          64 * KB_PER_GB, "SLURM_MEM_PER_CPU")
    # PER_CPU without a CPU count is not a budget; it must not be mistaken for one.
    check("SLURM_MEM_PER_CPU alone falls through",
          budget(tmp, environ={"SLURM_MEM_PER_CPU": "4096"}, meminfo_kb=500),
          500, "MemAvailable")
    # Explicit suffixes, which some sites set instead of a bare number.
    check("suffixed value",
          budget(tmp, environ={"SLURM_MEM_PER_NODE": "64G"}), 64 * KB_PER_GB,
          "SLURM_MEM_PER_NODE")
    # Garbage must not become a budget.
    check("unparseable value falls through",
          budget(tmp, environ={"PBS_RESC_MEM": "lots"}, meminfo_kb=777),
          777, "MemAvailable")


def test_cgroup_v2_leaf(tmp):
    root = os.path.join(tmp, "cg2leaf")
    write(os.path.join(root, "user.slice", "job-7", "memory.max"), str(48 * 1024**3))
    proc = os.path.join(tmp, "proc_cgroup_v2_leaf")
    write(proc, "0::/user.slice/job-7\n")
    check("cgroup v2 on the leaf",
          budget(tmp, cgroup_tree=root, proc_cgroup=proc), 48 * KB_PER_GB, "cgroup")


def test_cgroup_v2_ancestor(tmp):
    # The limit sits on the job scope, not on the leaf the process is in -- the case
    # a naive read of the mount root or of the leaf alone would miss entirely.
    root = os.path.join(tmp, "cg2anc")
    write(os.path.join(root, "pbs.slice", "job-7", "memory.max"), str(128 * 1024**3))
    write(os.path.join(root, "pbs.slice", "job-7", "task-1", "memory.max"), "max")
    proc = os.path.join(tmp, "proc_cgroup_v2_anc")
    write(proc, "0::/pbs.slice/job-7/task-1\n")
    check("cgroup v2 on an ancestor scope",
          budget(tmp, cgroup_tree=root, proc_cgroup=proc), 128 * KB_PER_GB, "cgroup")


def test_cgroup_v2_min_over_chain(tmp):
    # Two limits in the chain: the tighter one is the one that will actually kill us.
    root = os.path.join(tmp, "cg2min")
    write(os.path.join(root, "slice", "memory.max"), str(256 * 1024**3))
    write(os.path.join(root, "slice", "job", "memory.max"), str(64 * 1024**3))
    proc = os.path.join(tmp, "proc_cgroup_v2_min")
    write(proc, "0::/slice/job\n")
    check("cgroup v2 minimum over the chain",
          budget(tmp, cgroup_tree=root, proc_cgroup=proc), 64 * KB_PER_GB, "cgroup")


def test_cgroup_v2_unlimited(tmp):
    root = os.path.join(tmp, "cg2unl")
    write(os.path.join(root, "slice", "memory.max"), "max")
    proc = os.path.join(tmp, "proc_cgroup_v2_unl")
    write(proc, "0::/slice\n")
    check('cgroup v2 "max" is not a budget',
          budget(tmp, cgroup_tree=root, proc_cgroup=proc, meminfo_kb=4242),
          4242, "MemAvailable")


def test_cgroup_v1(tmp):
    root = os.path.join(tmp, "cg1")
    write(os.path.join(root, "memory", "pbspro", "job-9", "memory.limit_in_bytes"),
          str(96 * 1024**3))
    proc = os.path.join(tmp, "proc_cgroup_v1")
    write(proc, "8:memory,cpuacct:/pbspro/job-9\n2:cpu:/pbspro/job-9\n")
    check("cgroup v1 limit_in_bytes",
          budget(tmp, cgroup_tree=root, proc_cgroup=proc), 96 * KB_PER_GB, "cgroup")


def test_cgroup_v1_unlimited(tmp):
    # v1 writes PAGE_COUNTER_MAX rather than a sentinel string; treating that huge
    # number as a budget would be worse than falling through.
    root = os.path.join(tmp, "cg1unl")
    write(os.path.join(root, "memory", "memory.limit_in_bytes"),
          str(mod._CGROUP_UNLIMITED))
    proc = os.path.join(tmp, "proc_cgroup_v1_unl")
    write(proc, "8:memory:/\n")
    check("cgroup v1 PAGE_COUNTER_MAX is not a budget",
          budget(tmp, cgroup_tree=root, proc_cgroup=proc, meminfo_kb=1234),
          1234, "MemAvailable")


def test_meminfo_fallback(tmp):
    check("MemAvailable fallback", budget(tmp, meminfo_kb=8 * KB_PER_GB),
          8 * KB_PER_GB, "MemAvailable")


def test_nothing_readable(tmp):
    # The gate must degrade to "no budget" (pool sized by -c alone, as before the
    # gate existed) rather than to a bogus number.
    check("nothing readable", budget(tmp), None, "none")


def test_precedence(tmp):
    # A cgroup limit and a scheduler variable disagree: the scheduler wins, because
    # some sites enforce by polling and set no cgroup limit at all, while a cgroup
    # limit on a shared ancestor may be looser than the job's own allocation.
    root = os.path.join(tmp, "cgprec")
    write(os.path.join(root, "memory.max"), str(512 * 1024**3))
    proc = os.path.join(tmp, "proc_cgroup_prec")
    write(proc, "0::/\n")
    check("scheduler env beats cgroup",
          budget(tmp, environ={"SLURM_MEM_PER_NODE": "131072"},
                 cgroup_tree=root, proc_cgroup=proc, meminfo_kb=1600 * KB_PER_GB),
          128 * KB_PER_GB, "SLURM_MEM_PER_NODE")
    check("cgroup beats MemAvailable",
          budget(tmp, cgroup_tree=root, proc_cgroup=proc,
                 meminfo_kb=1600 * KB_PER_GB),
          512 * KB_PER_GB, "cgroup")


def capture_warning(source, environ):
    """Run the warning helper with stderr captured; return (fired, text)."""
    saved, buf = sys.stderr, io.StringIO()
    sys.stderr = buf
    try:
        fired = mod._warn_if_host_wide_budget(source, environ=environ)
    finally:
        sys.stderr = saved
    return fired, buf.getvalue()


def test_warning():
    fired, text = capture_warning("MemAvailable", {"PBS_JOBID": "1234.pbs"})
    assert fired, "warning did not fire on a host-wide budget inside a PBS job"
    assert "--max_memory" in text, "warning does not point at the flag: {!r}".format(text)
    assert "PBS_JOBID" in text, "warning does not name the scheduler: {!r}".format(text)
    print("  warning fires under a scheduler on a host-wide budget")

    for job_var in ("SLURM_JOB_ID", "SLURM_JOBID", "LSB_JOBID"):
        fired, _ = capture_warning("MemAvailable", {job_var: "1"})
        assert fired, "warning did not fire under {}".format(job_var)
    print("  warning fires for every recognised scheduler job variable")

    # Not under a scheduler: MemAvailable is a perfectly good budget, stay quiet.
    fired, text = capture_warning("MemAvailable", {})
    assert not fired and text == "", "warning fired on a plain host: {!r}".format(text)
    # Under a scheduler but with a real budget: nothing to warn about.
    for source in ("--max_memory", "AGENT_MEMORY", "PBS_RESC_MEM", "cgroup", "none"):
        fired, text = capture_warning(source, {"PBS_JOBID": "1234.pbs"})
        assert not fired and text == "", \
            "warning fired on source {}: {!r}".format(source, text)
    print("  warning stays quiet on a plain host and on a real budget")


def test_no_headroom_applied(tmp):
    """The helper returns the raw budget; the 0.8 headroom lives at the call site.

    TideCluster's equivalent applies its headroom inside the helper. If that were
    copied verbatim the factor would land twice (0.64) and quietly shrink the pool,
    so pin the raw value.
    """
    kb, _ = budget(tmp, explicit=100)
    assert kb == 100 * KB_PER_GB, "helper applied headroom: {} kB".format(kb)
    print("  helper returns the raw budget (no 0.8 headroom)")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        test_explicit_flag_wins(tmp)
        test_agent_memory(tmp)
        test_scheduler_vars(tmp)
        test_cgroup_v2_leaf(tmp)
        test_cgroup_v2_ancestor(tmp)
        test_cgroup_v2_min_over_chain(tmp)
        test_cgroup_v2_unlimited(tmp)
        test_cgroup_v1(tmp)
        test_cgroup_v1_unlimited(tmp)
        test_meminfo_fallback(tmp)
        test_nothing_readable(tmp)
        test_precedence(tmp)
        test_no_headroom_applied(tmp)
        test_warning()
    print("test_mem_budget: PASSED")
