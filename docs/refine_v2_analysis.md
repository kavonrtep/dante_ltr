# LTR boundary refinement — v2 design analysis

Status: empirical exploration; **no code changes yet**. This document records
the findings that should inform the next implementation of
`dante_ltr_refine`.

---

## 1. Background

The current `dante_ltr_refine` pipeline (added in commit `21d41f3`) performs
per-element LTR boundary refinement using parasail anchored extension across
clusters of similar 5'LTR sequences (mmseqs2 clusters within a lineage),
with optional MAFFT change-point fallback. In production we observed that
some refined boundaries lose the TG/CA dinucleotide motif and/or invalidate
the TSD that DANTE_LTR originally identified — i.e. the refinement worsens
the annotation in a meaningful fraction of cases.

This analysis quantifies the problem on two real datasets, tests an
alternative algorithm structure (use opposite-role LTRs as the anchor
pool), and reports the implications for a v2 design.

Datasets:

| dataset | source | LTR rows | clusters used (≥6 of each role) | rows in validation |
|---|---|---:|---:|---:|
| `at` | *Arabidopsis lyrata* (Alyr) | 5 174 | 85 | 2 826 |
| `g2` | *Pisum sativum*-like genome (g2) | 15 802 | 260 | 5 586 |

---

## 2. Definitions of terms

These are the conventions used throughout this document and proposed for
the v2 implementation.

**LTR roles.** Each complete element has two LTRs in biological orientation:

- **5'LTR** — the LTR at the 5' end of the element (in the element's own
  reading direction).
- **3'LTR** — the LTR at the 3' end of the element.

5'LTR and 3'LTR are *direct repeats*: in biological orientation both read
`TG…CA`.

**Boundaries (per LTR).** Each LTR has two boundaries:

- **Outer boundary** — the boundary between the LTR and adjacent genomic
  flank.
  - 5'LTR's outer boundary is biological 5' end → terminal `TG`.
  - 3'LTR's outer boundary is biological 3' end → terminal `CA`.
- **Inner boundary** — the boundary between the LTR and TE-internal
  sequence.
  - 5'LTR's inner boundary is biological 3' end → terminal `CA`.
  - 3'LTR's inner boundary is biological 5' end → terminal `TG`.

In both `dante_ltr_refine` v1 and the v2 design under discussion, **only
outer boundaries are refined** (the genomic-flank-adjacent ones).

**Motif.** The TG/CA dinucleotide present at every LTR's outer terminus
in biological orientation. `motif_ok = True` ↔ the 2 bp at the candidate
outer-boundary position matches `TG` (5'LTR) or `CA` (3'LTR).

**TSD (Target Site Duplication).** A short tandem duplication of genomic
sequence flanking the inserted element on both sides. A TSD looks like:

```
… genomic … TSD [TG-LTR-CA] [interior] [TG-LTR-CA] TSD … genomic …
```

DANTE_LTR detects TSD lengths of 4–8 bp (exact match preferred, ≤1
mismatch allowed at length ≥5). A TSD attribute `TSD=<seq>` is recorded
on every `transposable_element` row, and one `target_site_duplication`
GFF row is emitted on each side of the element.

**Cluster.** A group of LTRs from the same lineage clustered by mmseqs2
on 5'LTR sequence at a given identity (default 0.9, coverage 0.8). Sibling
3'LTRs are attached by parent.

**Anchor pools.** When refining the outer boundary of one LTR (the
*query*), the parasail anchored-extension algorithm aligns the query
against a set of *anchor* sequences. Anchors differ by which edge of
their LTR contributes to the alignment, and by what kind of flank that
edge has. The v2 terminology is:

- **outer-edge pool** (legacy name: A). Anchors are the OTHER same-role
  LTRs in the cluster. Each anchor contributes its OUTER edge — TG/CA
  with genomic flank — i.e. the same edge geometry as the query.
  *This is the current `dante_ltr_refine` (v1) behaviour.*
- **inner-edge pool** (legacy name: B). Anchors are the opposite-role
  LTRs in the cluster, used at their INNER edge — TG/CA with TE-internal
  flank. Same TG/CA terminus as the query (5'LTR and 3'LTR are direct
  repeats), but the TE-internal flank is *shared across cluster members*
  while the query's flank is *random genomic*; this asymmetry sharpens
  the boundary call.
- **joint pool** (legacy name: C). Union of outer-edge + inner-edge
  anchors. (P1 candidate from the design discussion.)

Where context is clear we write *outer / inner / joint*. The legacy
labels A / B / C are retained only in the validation TSV columns
(`A_g_snap`, `B_motif_ok`, …) because those files are already on disk.

| v2 name | legacy | role of anchors | anchor edge contributing | flank composition on anchor |
|---|---|---|---|---|
| outer pool | A | same-role  | outer | genomic |
| inner pool | B | opposite-role | inner | TE-internal |
| joint pool | C | both | outer + inner | mixed |

**Boundary estimate (`g`).** A 1-based genomic coordinate for the refined
outer boundary, output by parasail Nth-largest aggregation followed by
TG/CA snap within ±5 bp (the snap window).

**Inner-primary policy (formerly "v1 policy").** Rule for choosing the
boundary per side:

```
if  inner_motif_ok = TRUE  and  inner coordinate exists:
    use inner-pool coordinate
else:
    keep DANTE_LTR original coordinate
```

i.e. take the inner-pool call only when it passes the TG/CA motif check;
otherwise the side is unrefined. This is the policy whose simulations
appear in §4 (referred to there as "v1" because that is the variable
name used in the analysis scripts).

**Confidence labels (proposed).** A second pass labels each
inner-primary-accepted boundary by whether the outer-pool call agrees.
The label is informational; it does not change the boundary coordinate.

| label | condition |
|---|---|
| **dual**       | outer reached a coordinate, outer motif_ok, \|outer − inner\| ≤ 5 bp |
| **divergent**  | outer reached a coordinate, outer motif_ok, but \|outer − inner\| > 5 bp |
| **inner_only** | outer-pool call failed motif (or no call) |
| **unrefined**  | inner-pool also failed motif → DANTE_LTR original kept |

(Open to other keyword choices — see §6 Q7.)

**Pair-symmetry change Δ.** For an element with both LTRs:

```
W5 = width of 5'LTR (only the OUTER edge can shift in v1)
W3 = width of 3'LTR
asym_orig = |W5_orig − W3_orig|
asym_v1   = |W5_v1   − W3_v1|
Δ = asym_v1 − asym_orig
```

Δ > 0 means the v1 coords made the LTR pair more asymmetric than DANTE_LTR.

**TSD outcome.** At v1 (or refined) coordinates we re-detect the TSD using
the same algorithm DANTE_LTR uses (validated bit-for-bit against the
original on both datasets — see §4.4). Outcomes per element:

- **kept_exact** — same TSD sequence as the original GFF attribute.
- **kept_fuzzy** — original was 1-mismatch fuzzy; new still matches one of
  the original's two halves.
- **shifted** — both have a TSD but the sequence differs.
- **lost** — original had a TSD; v1 finds none.
- **gained** — original was `not_found`; v1 finds a TSD.
- **both_none** — neither side has a TSD.

---

## 3. What the validation script does

`tmp/validate_pool_refinement.py` (read-only) replays the existing
parasail anchored-extension geometry but runs it three times per
(member, side) — once with anchor pool A, once with B, once with C — and
records the resulting boundary coordinate, TG/CA status, snap offset,
and number of contributing alignments. Output: `tmp/pool_validation_at.tsv`
and `tmp/pool_validation_g2.tsv`.

`tmp/analyze_pair_and_tsd.py` re-detects TSDs at v1 boundaries (mirroring
the R algorithm), computes pair-symmetry deltas, and cross-tabulates
TSD outcome by which side(s) v1 changed.

The validation set is restricted to clusters with ≥6 of *both* roles (so
all three pools have enough voters); this is stricter than production
refine. Production data would have somewhat fewer eligible clusters; the
qualitative findings are robust to that.

---

## 4. Findings

### 4.1 Pool comparison — motif preservation and shift behavior

Original DANTE_LTR boundaries have TG/CA at 100.0% of rows in both datasets
(the input pipeline enforces this).

Column labels follow the validation TSV: `A` ≡ outer-edge pool,
`B` ≡ inner-edge pool, `C` ≡ joint pool.

| dataset | A motif | **B motif** | C motif | A=0 shift | B=0 shift | A vs B agree ≤5 bp |
|---|---:|---:|---:|---:|---:|---:|
| at | 67.4% | **92.8%** | 66.8% | 35.9% | 85.9% | 42.3% |
| g2 | 81.6% | **90.0%** | 82.3% | 57.0% | 80.8% | 71.0% |

**Key reading:** The inner-edge pool (B) preserves TG/CA on ~90 % of
rows and stays at the original coordinate on ~80–86 % of rows. The
outer-edge pool (A) loses motif on ~33 % of rows in `at` and ~18 % in
`g2`, and produces shifts > 20 bp on 28.6 % (at) / 9.1 % (g2) of rows.
The joint pool C tracks A almost exactly — Nth-largest aggregation lets
the longer same-role extensions dominate.

**Why the inner pool works.** The inner-edge anchor's flank is
TE-internal sequence *shared* across cluster members; the query's flank
is *random genomic*. Parasail extension cannot extend conservation past
the LTR terminus because the flanks don't match. The outer pool's
anchor flanks are both random genomic — chance flank similarity drags
the consensus extension past TG in ~30 % of cases.

**Dominant regression pattern** (orig motif `T`, outer pool `F`, inner pool `T`):

| | at | g2 |
|---|---:|---:|
| pattern count | 770 (27.2 %) | 521 (9.3 %) |
| total motif loss by outer pool | 921 (32.6 %) | 1 027 (18.4 %) |
| of those, inner pool carries motif | 785 (85.2 %) | 855 (83.3 %) |

When the outer pool loses motif, the inner pool preserves it ~84–85% of
the time on the same row.

### 4.2 Inner-primary policy simulation (revert if no motif)

| dataset | final motif | rows changed | shift > 5 bp |
|---|---:|---:|---:|
| at | **100.0 %** | 9.5 % | 6.7 % |
| g2 | **100.0 %** | 11.6 % | 8.6 % |

Refinement fires on ~10 % of rows; the remaining ~90 % keep DANTE_LTR's
coordinate. No motif is ever lost (by construction).

Confidence-label distribution for all rows (using the new keywords):

| label | at | g2 |
|---|---:|---:|
| **dual** (outer ∩ inner agree, both motif) | 39.5 % | 63.8 % |
| **divergent** (outer motif, disagrees > 5 bp) | 25.5 % | 15.8 % |
| **inner_only** (outer motif fail) | 27.8 % | 10.4 % |
| **unrefined** (inner no motif → kept original) | 7.2 % | 10.0 % |

For the ~10 % of rows where the inner-primary policy actually fires:

| breakdown | at (n=268) | g2 (n=648) |
|---|---:|---:|
| outer also agrees ≤ 5 bp | 43.3 % | 65.6 % |
| outer motif but disagrees > 5 bp | 34.0 % | 22.8 % |
| outer motif fail | 22.8 % | 11.6 % |

### 4.3 Pair symmetry under inner-primary policy (per-element)

| dataset | elements with both LTRs | same | grew | shrank | grew > 5 bp |
|---|---:|---:|---:|---:|---:|
| at | 1 413 | 82.1 % | 17.1 % | 0.8 % | **12.2 %** |
| g2 | 2 793 | 79.0 % | 20.3 % | 0.7 % | **14.1 %** |

**Inner-primary *grows* asymmetry on ~17–20 % of pairs**; on ~12–14 %
the increase is more than 5 bp. This is because the policy evaluates
each side independently — when only one side gets an inner-pool
correction, the unmoved side stays at DANTE_LTR's coord and the pair
loses its near-equal-width property.

### 4.4 TSD detection — sanity check against R

Re-detecting TSD at *original* coords using the Python port (in
`analyze_pair_and_tsd.py`) and comparing to the input-GFF `TSD=` attribute:

| dataset | seq match | both none | only-GFF | only-Py | seq differ |
|---|---:|---:|---:|---:|---:|
| at | 72.6 % | 27.4 % | 0 | 0 | 0 |
| g2 | 68.7 % | 31.1 % | 0 | 0.1 % (4 elements) | 0 |

The Python implementation reproduces R's TSD detection bit-for-bit
modulo a handful of fuzzy edge cases. We can use the Python detector for
v2 implementation without an R subprocess.

### 4.5 TSD outcome at inner-primary boundaries

For all elements with both LTRs in the validation set:

| outcome | at | g2 |
|---|---:|---:|
| kept_exact | 60.7 % | 59.5 % |
| kept_fuzzy | 4.5 % | 2.6 % |
| shifted | 1.4 % | 0.9 % |
| **lost** | **5.9 %** | **5.8 %** |
| gained | 4.5 % | 5.9 % |
| both_none | 22.9 % | 25.3 % |

Restricted to elements where the inner-primary policy actually changed at least one side:

| outcome | at (n=253) | g2 (n=591) |
|---|---:|---:|
| kept_exact | 1.6 % | 0.3 % |
| **lost** | **33.2 %** | **27.4 %** |
| gained | 25.3 % | 27.4 % |
| both_none | 31.6 % | 40.6 % |
| shifted | 7.9 % | 4.2 % |

**When the inner-primary policy fires, the existing TSD is destroyed
~30 % of the time.** This matches the failure mode the user observed
in production.

### 4.6 The decisive cross-tab — coupling × TSD

Splitting elements by *which* side the inner-primary policy changed:

| change pattern | at: count, % | TSD kept | TSD lost | TSD gained | TSD both_none |
|---|---:|---:|---:|---:|---:|
| no change | 1 160 (82.1 %) | 79.1 % | — | — | 20.9 % |
| only 5'LTR | 162 (11.5 %) | 0.6 % | **32.1 %** | 26.5 % | 34.6 % |
| only 3'LTR | 76 (5.4 %) | 5.3 % | **38.2 %** | 26.3 % | 17.1 % |
| both sides | 15 (1.1 %) | 0 % | 20.0 % | 6.7 % | 73.3 % |

| change pattern | g2: count, % | TSD kept | TSD lost | TSD gained | TSD both_none |
|---|---:|---:|---:|---:|---:|
| no change | 2 202 (78.8 %) | 78.6 % | — | 0.2 % | 21.2 % |
| only 5'LTR | 299 (10.7 %) | 0.0 % | **34.4 %** | 23.7 % | 37.8 % |
| only 3'LTR | 235 (8.4 %) | 0.9 % | **23.0 %** | 36.2 % | 34.5 % |
| both sides | 57 (2.0 %) | 0 % | 8.8 % | 10.5 % | 80.7 % |

**One-sided inner-primary changes destroy the original TSD in
23–38 % of the cases where one existed.** No-change pairs preserve
TSD ~79 %. Both-sided changes are rare (1–2 %) and almost never
re-create the original TSD, suggesting they reach into territory
DANTE_LTR didn't validate.

---

## 5. What the data says about v2 design

Three things are now empirically supported:

1. **Use the inner-edge pool as the boundary signal**, not the outer-edge
   pool. The inner pool's flank-asymmetry property gives a sharp
   LTR-vs-non-LTR boundary call; the outer pool drifts. The joint pool
   does not work as originally hoped — it inherits the outer-pool drift.

2. **The outer-edge pool is a validator, not a primary.** It produces
   clean confidence labels (§4.2) on the inner-pool calls. The
   `dual`-confidence bucket is the modal outcome on real data
   (40 % at, 64 % g2).

3. **Per-side acceptance is not enough — pair coupling matters.** Even
   when the inner-primary policy preserves motif on each side
   independently (100 %), one-sided changes destroy TSD on ~30 % of
   those events (§4.6). This is exactly the failure mode the user
   reported in production.

The cleanest rule that unifies all three signals:

> **Accept an inner-primary change on a side only if (a) the inner
> pool's motif is OK on that side AND (b) accepting it does not
> invalidate the original TSD, OR a coupled change on both sides
> preserves/recreates a TSD.**

A coupled change is allowed when both sides have inner-pool proposals
AND the new flanks support a TSD (kept_exact, kept_fuzzy, shifted, or
gained). One-sided changes that destroy the original TSD revert.

---

## 6. Design decisions

Decisions captured in this round (after the analysis, before the
implementation plan):

1. **Pair-coupling level — per-side or all-or-nothing per element?**
   ✅ **Per-side accept/revert (Option α).** Each side is gated
   independently on motif + per-element TSD recheck. We do not
   require both sides to refine together; we require that whatever
   refinement is accepted preserves the TSD context.

2. **TSD-loss revert: AND-rule or OR-rule with motif?**
   ✅ **TSD-only test (collapsed AND-rule).** Under the inner-primary
   policy the motif is always preserved on accepted sides by
   construction, so `motif_lost AND TSD_lost` collapses to
   `TSD_lost`. The rule is: *if applying an inner-pool refinement to
   a side would destroy the originally-detected TSD, revert that
   side.*

3. **Elements where DANTE_LTR found no original TSD (`both_none`).**
   ✅ **TSD-loss rule does not apply** (no TSD to defend; only the motif
   test gates the change, which the inner-primary policy already
   satisfies). A *gained* TSD on a one-sided change is treated as a
   weak positive confirmation of the refinement, not a spurious
   match.

4. **MAFFT fallback.**
   ✅ **Keep as a safety net, but adapt to the inner-primary logic.**
   Cluster validation rate should be computed against the inner pool
   (which has the higher pass rate), and the fallback should fire
   only on clusters where the inner pool fails — not where the outer
   pool fails (which is now expected and not a problem). Detailed
   wiring TBD in the implementation plan.

5. **Joint-MSA refinement (option C-3 from the v2 plan).**
   ✅ **Out of scope for v2; defer to v3.** Keep the design idea
   documented for later. v2 sticks to per-element parasail with the
   inner-edge pool as the primary signal.

6. **Cluster-size threshold.**
   ✅ **Require ≥ N of each role** (currently default `N = 6`). This
   matches the validation script's filter and ensures both pools
   have enough voters.

7. **Confidence labels emitted as GFF attribute?**
   ✅ **Yes, emit; rename to be more informative.** Initial proposal:
   `dual / divergent / inner_only / unrefined` (defined in §2). The
   GFF attribute name will be `Refinement_Confidence` (carried over
   from v1). Naming is open to revision in the implementation plan.

8. **Pool naming — replace A / B / C.**
   ✅ **Renamed to `outer-edge / inner-edge / joint pool`** (or
   *outer / inner / joint* in shorthand). This is what the v2 code
   and docs use; the A/B/C labels survive only in the validation
   TSV column headers already on disk.

9. **Confidence-label keywords.**
   ✅ **`dual / divergent / inner_only / unrefined`** (defined in §2).

10. **TSD-gain attribution in GFF3.**
    ✅ **Update both the TE `TSD=` attribute and the per-side
    `target_site_duplication` child rows.** When an inner-primary
    change gains a new TSD, replace the input `TSD=not_found` with
    the new TSD sequence, and emit two new `target_site_duplication`
    rows as children of the TE. Child rows inherit the parent's
    `Refinement_Confidence` label.

---

## 7. Open follow-ups

These are still open and need a decision after v2 ships:

- **MAFFT fallback wiring.** What "inner-pool failure rate" threshold
  triggers it; whether fallback proposals also have to pass the
  TSD-loss revert rule. *Deferred — we may do some testing on v2
  output before re-tuning the fallback.*

---

## 8. Reproducibility

```bash
conda activate dante_ltr

python tmp/validate_pool_refinement.py \
    -g test_data/at/Alyr_dante_ltr.gff3 \
    -s test_data/at/Alyr.fasta \
    -o tmp/pool_validation_at.tsv \
    --threads 6

python tmp/validate_pool_refinement.py \
    -g test_data/g2/DANTE_LTR.gff3 \
    -s test_data/g2/genome.fasta \
    -o tmp/pool_validation_g2.tsv \
    --threads 6

python tmp/analyze_pair_and_tsd.py
```

Outputs:

- `tmp/pool_validation_at.tsv`, `tmp/pool_validation_g2.tsv` — per-row
  estimates A/B/C with motif and snap diagnostics.
- `tmp/pool_validation_at_run.log`, `tmp/pool_validation_g2_run.log` —
  validation run summary printed to stderr.

The two analysis scripts (`tmp/validate_pool_refinement.py`,
`tmp/analyze_pair_and_tsd.py`) are read-only and **not** intended for
inclusion in the production code base.
