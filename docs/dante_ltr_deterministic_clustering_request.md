# Feature request (dante_ltr): deterministic library clustering

**Status:** resolved in two parts — canonical input sort, then `--spaced-kmer-mode 0` (see Resolution).
**Affects:** `dante_ltr_to_library` → `utils/mmseq_clustering.R` (and, second-order, `utils/extract_fasta.R`).
**CARP tracking:** part of `fix/deterministic-library-clustering`.

## Summary

`dante_ltr_to_library` clusters the extracted transposable-element sequences
(`TE_all.fasta`) with `mmseqs easy-cluster`, and takes the cluster
representatives verbatim as the LTR-RT library. **`mmseqs easy-cluster` is
order-sensitive**: the same set of input sequences presented in a different
order yields a *different set of representative consensi and a different cluster
count*. Because `TE_all.fasta` is emitted in the row order of `DANTE_LTR.gff3`
(rtracklayer `import` preserves file order), and that order is not canonical
(see below), the LTR-RT library — and therefore the downstream RepeatMasker
annotation — is not reproducible run-to-run / across machines.

This is the same class of bug, with the same fix, that TideCluster already
addressed in its comparative analysis (TideCluster `changelog.md`, 1.13.1,
issue #4: "Deterministic comparative analysis").

## Evidence (measured in CARP)

- `mmseqs easy-cluster` on a fixed LINE input, same file, at 1 and 16 threads →
  **byte-identical** representatives. So it is **not** thread-related.
- The *same* sequences **shuffled** → **~19 % of representatives differ** and the
  **cluster count changes** (e.g. 2250 → 2260). Re-sorting the input into one
  canonical order before clustering restores byte-identical representatives.
- On one genome (GCA_964200825.2) the Ty1_copia/Angela library consensi differ
  ~49 % between pipeline runs; the elements themselves are stable (same count),
  so the churn is representative-election over a stable set — i.e. exactly this
  order-sensitivity.
- `DANTE_LTR.gff3` row order is **not canonical** on the multi-chunk path: chunk
  results are concatenated in chunk-index order, and the chunk count depends on
  genome size / open-file limit / machine, so `TE_all.fasta` order (hence the
  clustered library) varies across environments.

## Requested fix (in dante_ltr)

Make the library a deterministic function of the input **set**, independent of
record order. Preferred: sort `TE_all.fasta` into a canonical order **by
sequence content** immediately before clustering, inside
`utils/mmseq_clustering.R` (so every caller of `dante_ltr_to_library` benefits).
An out-of-core sort keeps it cheap on large inputs, e.g.

```sh
seqkit fx2tab TE_all.fasta \
  | LC_ALL=C sort -t$'\t' -k2,2 -k1,1 -S <buf> --parallel=<threads> -T <tmp> \
  | seqkit tab2fx > TE_all.sorted.fasta
```

Sorting by **sequence** (not coordinate/ID) makes the order invariant to
upstream coordinate jitter and to chunk grouping. Optionally also expose a
`--deterministic` switch that additionally forces `mmseqs --threads 1` for a
byte-identical `*_rep_seq.fasta` (matching TideCluster's flag), though sorting
alone already fixes the representative set.

Performance note: the sort is O(n log n) over bytes the clustering reads anyway
and is disk-backed, so it is a small fraction of clustering wall-time even on
30–90 Gbp genomes with large TE sets.

## What CARP does in the meantime

CARP cannot fix this without modifying the dependency, so it does **not** work
around it inside `make_library_of_ltrs`. CARP *does* independently sort the
inputs to the clustering steps it drives itself (`dante_line`, `reduce_library`
CAP3/mmseqs, `make_tir_combined_library`), so the LTR consensi are re-clustered
deterministically at the `reduce_library` stage. Full reproducibility of the
Angela/LTR layer requires this dante_ltr fix to land as well.


## Resolution

Fixed in two independent parts, both in `utils/mmseq_clustering.R`.

**Part 1 — input order.** `TE_all.fasta` is sorted canonically by sequence
content (then name) with `order(..., method = "radix")` immediately before
partitioning, which is the fix requested above. The in-memory sort replaced the
proposed `seqkit`/`sort` pipeline: the sequences are already read into memory to
be partitioned, so the out-of-core route added a dependency and a temp file for
no benefit.

**Part 2 — random spaced k-mer pattern.** The sort alone was *not* sufficient.
Measured afterwards: with a byte-identical `partitioned_s900_w1000.fasta`,
repeated `mmseqs easy-cluster` runs still produced 5146 / 5162 / 5167 / 5169 /
5173 clusters. Bisected to `--spaced-kmer-mode 1`, which linclust passes
explicitly; with no `--spaced-kmer-pattern` supplied mmseqs generates a random
pattern per process. The first divergent intermediate is `linclust/pref`,
kmermatcher's output, and everything downstream inherits it. Passing
`--spaced-kmer-mode 0` makes it deterministic.

Note for the evidence above: the observation that 1 and 16 threads gave
byte-identical representatives on the LINE input was luck, not a property —
this failure mode varies at `--threads 1` too, and is unaffected by
`--split-memory-limit`. The `--deterministic` switch forcing `--threads 1`
would therefore *not* have fixed it; it is not needed and was not added, since
the output is now identical across thread counts anyway.
