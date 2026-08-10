# Implementation plan: parallelize the chunk loop in `dante_ltr`

Companion to `docs/dante_ltr_chunk_parallelism_request.md`. The request is valid; this
plan turns the serial per-chunk loop (`dante_ltr:1246–1257`) into a memory-gated process
pool of single-threaded children, with byte-identical output.

## Goals / non-goals

- **Goal:** run the independent chunk subprocesses concurrently, bounded by cores *and*
  available memory, without changing output.
- **Goal:** no orphaned R processes on failure or Ctrl-C.
- **Non-goal (this pass):** the whole-genome else-branch (`dante_ltr:1312`) is untouched —
  it legitimately keeps `-c args.cpu`.
- **Non-goal (this pass):** fixing the pre-existing run-to-run RNG nondeterminism
  (`ltr_utils.R:534` `sample()` under `mclapply` with no global seed) and per-chunk
  `TE_########` ID collisions. Both are orthogonal and unchanged by this work.

## Why output stays identical

- Concatenation (`:1299`), GFF back-remap (`:1280`), and stats sum all iterate
  `range(number_of_temp_files)` in **index order**, not completion order — so
  `imap_unordered` cannot reorder anything.
- `get_unique_features` (`:621`) re-sorts by `(start, header)`, independent of line order.
- Children switch from `-c args.cpu` to `-c 1`. `mclapply` uses
  `mc.set.seed=TRUE, mc.preschedule=FALSE`, which seeds per element index independent of
  `mc.cores`, so per-chunk output is expected core-count-invariant. **Verify** with the
  determinism diff in the Testing section rather than assume.

## Design

### 1. Module-level worker (picklable)

`multiprocessing.Pool` requires a top-level function. Add near the other helpers:

```python
def _run_detect_ltr(cmd):
    """Run one detect_putative_ltr.R invocation. Raises CalledProcessError on failure
    (picklable, so it propagates out of the pool)."""
    subprocess.check_call(cmd)
```

### 2. Factor out the per-chunk command

Build the command in one place so the probe run and the pooled runs are identical:

```python
def _detect_ltr_cmd(tool_path, fasta, gff, out, cpu, args, constrains_path):
    cmd = [f'{tool_path}/utils/detect_putative_ltr.R',
           '-s', fasta, '-g', gff, '-o', out,
           '-c', str(cpu), '-M', str(args.max_missing_domains),
           '-L', str(args.min_relative_length)]
    if constrains_path:
        cmd += ['--te_constrains', constrains_path]
    if args.no_ambiguous_domains:
        cmd += ['--no_ambiguous_domains']
    return cmd
```

### 3. Longest-first ordering

Chunk files are already balanced by round-robin over size-sorted ids, but a few huge
scaffolds can unbalance them. Order chunk indices by input FASTA size descending. This
also gives us the **largest chunk as the memory probe** for free:

```python
order = sorted(range(number_of_temp_files),
               key=lambda i: os.path.getsize(temp_files_fasta[i]), reverse=True)
```

### 4. Memory-gated pool size (measure at `-c 1`, TideCluster pattern)

Do **not** anchor on the request's 108.9 GB — that was measured with `-c 96` children,
whose `mclapply` fork fan-out inflates RSS. Measure the real single-threaded footprint:

```python
def _mem_available_kb():
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1])   # kB
    except OSError:
        pass
    return None

# Probe: run the largest chunk serially with -c 1, then read its peak RSS.
probe = order[0]
subprocess.check_call(_detect_ltr_cmd(tool_path, temp_files_fasta[probe],
                                      temp_files_gff[probe], output_files[probe],
                                      1, args, constrains_path))
peak_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss  # largest child, kB on Linux

remaining = order[1:]
avail = _mem_available_kb()
mem_cap = args.cpu
if avail and peak_kb > 0:
    mem_cap = max(1, int(avail * 0.8) // peak_kb)   # 0.8 safety margin
pool_size = max(1, min(args.cpu, mem_cap, len(remaining)))
if pool_size < min(args.cpu, len(remaining) + 1):
    print('memory-gating chunk pool to {} concurrent chunks '
          '(~{} MB/chunk, {} MB available)'.format(
              pool_size, peak_kb // 1024, (avail or 0) // 1024))
```

`RUSAGE_CHILDREN.ru_maxrss` after the probe reflects that one child's peak (largest child,
not a sum). `getpagesize`/units are kB on Linux — the target platform.

### 5. Run the remaining chunks in the pool

```python
from multiprocessing import Pool

cmds = [_detect_ltr_cmd(tool_path, temp_files_fasta[i], temp_files_gff[i],
                        output_files[i], 1, args, constrains_path) for i in remaining]
if cmds:
    with Pool(pool_size) as p:
        for _ in p.imap_unordered(_run_detect_ltr, cmds):
            pass
```

`with Pool(...)` calls `terminate()` on exception, so a `CalledProcessError` from any
child tears the pool down instead of leaking R processes. Add `import` at top of file
(`from multiprocessing import Pool`); `resource` is already imported.

### 6. File-descriptor headroom must scale with pool width

`_limit_chunks_to_fd_budget(number_of_temp_files, fd_headroom=64)` currently assumes a
single child. With `pool_size` concurrent R processes each opening their own handles, bump
the reserve. Since pool size is only known after chunking, keep the split-time budget as
is but pass a headroom that accounts for concurrency, or simply raise the default reserve
to cover `args.cpu`:

```python
number_of_temp_files = _limit_chunks_to_fd_budget(
    number_of_temp_files, fd_headroom=64 + args.cpu)
```

(Each concurrent child needs a handful of descriptors for its temp files/BLAST; `+args.cpu`
is a cheap over-estimate.)

### 7. Signal handling (hardening)

To avoid orphaned R children on Ctrl-C, have pool workers ignore `SIGINT` and let the main
process handle it and call `terminate()`:

```python
import signal
def _pool_worker_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
# Pool(pool_size, initializer=_pool_worker_init)
```

Wrap the pool block in `try/except KeyboardInterrupt: p.terminate(); raise`.

## Optional phase 2: reclaim idle cores when chunk-count-bound

When `pool_size` is limited by **chunk count** (`len(remaining) < args.cpu`) rather than
by memory, cores sit idle at `-c 1`. Only then, give each child extra threads:

```python
threads = 1
if pool_size == len(remaining) and pool_size < args.cpu:
    threads = max(1, args.cpu // pool_size)
```

Do **not** do this when memory-bound — extra threads mean `mclapply` forking, which raises
per-child RSS and breaks the gate. Skip phase 2 entirely if it complicates review; the
massively-chunked large-genome case (the one the request is about) is memory/chunk-rich and
gains nothing from it.

## Docs / help text

- Update `-S/--max_chunk_size` help: it now controls per-chunk memory **and** pool
  granularity (chunk count = pool work items).
- Note in README / CLAUDE.md pipeline section that per-chunk detection now runs in a
  memory-gated pool. Keep it to one or two sentences (repo doc-style convention).

## Testing

1. **Determinism (the critical check).** On `test_data`, force chunking with a small `-S`
   so `total_size > max_chunk_size and number_of_sequences > 1`, run the serial `main`
   (current code) and the pooled `main`, and `diff` `.gff3`, `_statistics.csv`, and each
   `_D*/DL*` FASTA. Must be byte-identical. Also run the pooled path twice to see the
   pre-existing RNG variance (expected, orthogonal — document, don't fix here).
2. **Core-invariance of a chunk.** Run one chunk directly via `detect_putative_ltr.R` at
   `-c 1` vs `-c 8`; diff outputs to confirm the `-c 1` child switch is safe.
3. **Memory gate.** Log the probe `peak_kb` and chosen `pool_size` on a chunked run;
   sanity-check `pool_size × peak_kb < MemAvailable`.
4. **Failure propagation.** Point one chunk at a malformed GFF; confirm the run aborts with
   the child's error and no R processes remain (`pgrep -f detect_putative_ltr`).
5. `./tests.sh` end-to-end.

## Rollout

- Small, self-contained change confined to the chunked branch of `main()` plus two helpers.
- Bump `version.py` per repo convention.
- Reproduce the memory-gate math locally on a scratch multi-sequence FASTA before relying
  on CI (per the "reproduce slow CI locally" convention).
