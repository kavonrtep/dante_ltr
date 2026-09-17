# `dante_ltr_to_library`: core mode's mixed-precision labels drop 15% of clusters and discard recoverable lineage calls

DANTE_LTR 0.6.0.0's `--mode core` is exactly what CARP needs for genomes far from
REXdb — on a *Draparnaldia* assembly it turns 0 complete elements into 1,408. CARP
plans to expose it as an opt-in `dante_ltr_mode: core`
(see [`dante_ltr_core_mode_plan.md`](dante_ltr_core_mode_plan.md)).

Building the repeat library from that output surfaced one interaction that was
invisible in lineage mode: **`dante_ltr_to_library` treats "labelled at a coarser
depth" as a classification conflict and discards the cluster.** Core mode produces
such labels by design, so the rule fires constantly — and in the process it also
throws away lineage evidence that only becomes visible at cluster level.

Observed with `dante_ltr` 0.6.0.0 and the bundled `utils/mmseq_clustering.R`,
default parameters (`-m 3 -p 0.95`).

## Finding 1 — an ancestor label is treated as a conflict

`utils/mmseq_clustering.R` decides each cluster's fate with:

```r
final_name      <- ifelse(main_class_proportion > opt$proportion_min,
                          main_class_name, consensus_names)
resoved_names_l <- final_name == main_class_name
final_name      <- final_name[resoved_names_l]        # everything else is dropped
```

When no label holds >95% of the members, `resolve_name()` computes the LCA of the
distinct labels, and the cluster survives **only if that LCA happens to equal the
majority label**. A cluster whose members are `Class_I|LTR|Ty3/gypsy` and
`Class_I|LTR|Ty3/gypsy|chromovirus` has LCA `Ty3/gypsy` and majority
`chromovirus`, so it is discarded — although the two labels do not contradict each
other at all. One is simply less resolved than the other, which is the normal
output of core mode's LCA-based classification.

Measured on the *Draparnaldia* core-mode run, re-deriving every cluster decision
from `mmseqs_cluster.tsv`:

```
clusters with >=3 members            464
  kept, >95% majority                355
  kept, LCA == majority               38
  DROPPED by the rule above           71
     distinct labels form a single
     ancestor chain (no conflict)     66   <- 93% of the loss
     sibling lineages                  0
     cross-superfamily (copia/gypsy)   5   <- genuine conflicts
```

The same re-derivation on this genome's **lineage-mode** output: 13 clusters, 0
dropped. The rule is invisible until core mode produces internal-node labels.

## Finding 2 — cluster-level lineage evidence is discarded

Core mode assigns an element the LCA of *its own* domains' calls, so an element
whose RT matches Tekay while RH/INT/CHD match Chlamyvir is labelled
`chromovirus`. That is correct per element. But a cluster pools ~40 elements, and
a cluster may contain 18 elements that each independently resolved to Chlamyvir
and none that resolved to any other lineage. That is evidence no single element
had — the cluster is the family, so it is the right level at which to make the
call. Today the cluster is labelled with the LCA and the lineage information is
dropped.

Measured, over clusters that survive (or would survive under Request 2), whose
distinct labels form a single chain with exactly one deepest label:

```
deepest label carried by >=2 distinct elements                  88 clusters
  ... and by >=25% of the cluster's elements                    69 clusters
        promoted to Class_I|LTR|Ty3/gypsy|chromovirus           50
        promoted to Class_I|LTR|Ty3/gypsy|chromovirus|Chlamyvir 19
```

Choosing the share threshold — distribution of (elements carrying the deepest
label) / (elements in the cluster), over those same 88 clusters:

```
quantiles   10% = 0.15   25% = 0.31   50% = 0.67   75% = 0.79
>= 0.10 -> 86    >= 0.20 -> 71    >= 0.25 -> 69    >= 0.50 -> 62
```

The choice is not sensitive: 0.25 and 0.50 differ by 7 of 88 clusters. 0.25 is
proposed because the 25th percentile sits at 0.31, so it trims the tail of
single-stray-element cases without cutting into the bulk.

Note the *other* candidate denominator — share of the already-resolved members —
is degenerate: it is 1.00 from the 5th percentile up, because almost every chain
is two-level, so "resolved" and "deepest" are the same set. Share of **all**
members is the ratio that discriminates.

## Impact on downstream consumers

For CARP the cluster representatives *are* the LTR library handed to
RepeatMasker, so a dropped cluster is unmasked genome.

```
library on the Draparnaldia core-mode run, as written today
  397 sequences / 435,672 bp
re-derived with Requests 2+3
  459 sequences / 503,834 bp      (+66 clusters, +17% bp)
```

Nineteen of those clusters would additionally carry a real `Chlamyvir` call
instead of the bare `chromovirus` bucket, which flows straight through to the
RepeatMasker annotation, `summary_statistics.csv` and the per-class density
tracks.

## Requests

**R1 — a flag, defaulting to today's behaviour.**
`dante_ltr_to_library --annotation_conflict {strict,nested}`, default `strict`,
byte-identical to the current code path. CARP would pass `nested` only when
running `--mode core`. It must stay opt-in rather than keyed off the mode:
lineage mode can also emit internal-node element labels (e.g.
`Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Tat`), so `nested` is not provably a
no-op there on every genome.

**R2 — under `nested`, recover ancestor/descendant chains.**
Keep the cluster when its distinct labels form a single ancestor chain — every
pair ancestor-or-descendant — and label it with the shallowest, which is the LCA
and is itself one of the labels present. Sibling-lineage and cross-superfamily
mixes stay dropped, exactly as today. The stricter single-chain test is preferred
over "the LCA is present" because it can never place two named lineages under one
representative; on this genome the two formulations are identical.

**R3 — under `nested`, count distinct source elements, not sliding windows.**
Members of a cluster are 1 kb windows, so one long element contributes ~8 members
and outvotes several short ones — which matters for R4's share threshold. Collapse
`_sliding:` suffixes to the element id before counting.
**This must apply only under `nested`:** re-deriving `strict` with element
collapse changes its outcome on 4 of 464 clusters here (393 vs 397 kept), so
`strict` has to keep counting windows to stay byte-identical.

**R4 — under `nested`, promote to the deepest label when the evidence supports it.**
After R2, if exactly one deepest label exists, it is carried by at least
`--lineage_promotion_min_elements` (default 2) distinct elements, and those are at
least `--lineage_promotion_min_share` (default 0.25) of the cluster's elements,
label the cluster with that deepest label instead of the LCA. Promotion belongs
inside `nested` rather than behind its own flag: 27 of the qualifying clusters are
ones `strict` already keeps, so an independent switch would change default output.

Two invariants worth asserting in code, because they are what makes this safe:
the promoted label is always a descendant of the LCA (refinement along one path,
never a reclassification), and promotion never fires when two sibling lineages are
present.

**R5 — log one summary line** — clusters kept / recovered / promoted / dropped —
so the effect is visible without re-deriving it from `mmseqs_cluster.tsv`.

## Reproduction

```bash
# 1. core-mode detection
dante_ltr --mode core -g DANTE.gff3 -s genome.fasta -o core -c 12 --max_memory 48

# 2. library build, default parameters
dante_ltr_to_library --gff core.gff3 --output_dir lib_core -s genome.fasta -c 8

# 3. the decisions are all re-derivable from
#    lib_core/mmseqs2/mmseqs_cluster.tsv   (member -> cluster representative)
#    lib_core/mmseqs2/mmseqs_rep_seq.fasta (representative sequences)
#    Each member name is  <element_id>#<Final_Classification>_sliding:<start>-<end>
```

The same two commands with the default `--mode lineage` give 13 clusters, none
dropped — the contrast that localises the finding to core mode's label depth
rather than to the clustering itself.

## Evidence from the CARP side

`Class_I|LTR|Ty3/gypsy`, `Class_I|LTR|Ty3/gypsy|chromovirus` and
`Class_I|LTR|Ty3/gypsy|non-chromovirus` are already in CARP's
`classification_vocabulary.yaml` as valid internal nodes, and all 25 distinct
classifications core mode emitted across these runs pass CARP's
`classification.py validate`. Nothing about the coarse labels is malformed — they
are a legitimate, deliberate output of core mode, and only the library builder's
conflict rule treats them as a problem.
