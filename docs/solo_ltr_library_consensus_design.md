# Improved LTR Consensus Generation for `dante_ltr_solo`

**Status:** design draft, no code changes yet
**Scope:** `utils/build_ltr_library.R` (the per-cluster consensus step) and the clustering logic that feeds it
**Goal:** generate LTR library consensus sequences whose boundaries reflect the *true* LTR/internal junction, by leveraging the asymmetric conservation pattern of flanking sequences in pooled 5'LTR + 3'LTR alignments.

---

## 1. Background — current behaviour

Today (`utils/build_ltr_library.R:142-238`), the library is built like this:

1. All `long_terminal_repeat` features are extracted from the DANTE_LTR GFF3 — both 5'LTRs and 3'LTRs together — using **exact feature coordinates with no flanking** (`getSeq(s, ltr_feat)` at line 164). The `--flank` option is declared but never used.
2. LTRs are partitioned by `Final_Classification` (lineage).
3. Within each lineage, MMseqs2 (`easy-cluster`, `--min-seq-id 0.9 -c 0.8`) groups LTRs into clusters. There is **no constraint that the 5'LTR and 3'LTR of the same element end up in the same cluster.**
4. For each cluster of size ≥ 3, MAFFT `--auto` produces an MSA, and a column-wise majority consensus is built. Columns with > 50% gaps are dropped (this is the only "trimming"). Clusters smaller than 3 fall back to the highest-rank single member, no alignment.
5. The resulting consensus has no concept of "LTR boundary" beyond what the input annotation already declared — if the original DANTE_LTR boundary was wrong by N bp, the consensus is wrong by N bp too.

### What this misses

- **No use of flanks.** Boundary information that lives just outside the annotated LTR is thrown away.
- **No element-pair coupling.** A 5'LTR and its sibling 3'LTR can be split across MMseqs2 clusters when one of them is degraded, even though biologically they are the same sequence at insertion time.
- **No boundary correction.** The consensus inherits whatever errors the input annotation has. Subsequent BLAST search in `detect_solo_ltr.R` then propagates those errors into every solo-LTR call.
- **No cross-validation between 5' and 3' LTRs.** The two ends of an element are near-identical at insertion, so they should produce mutually consistent boundaries — that consistency is currently unused.

---

## 2. Key insight — asymmetric flank conservation

Once we extract LTRs **with flanking sequence**, the flanks contain a strong, asymmetric, exploitable signal:

| LTR copy   | 5' flank (just upstream)                    | LTR body              | 3' flank (just downstream)                |
|------------|---------------------------------------------|-----------------------|-------------------------------------------|
| **5' LTR** | random genomic — *non-conserved across copies* | conserved (LTR family) | TE internal start (PBS, gag) — *conserved across copies of same family* |
| **3' LTR** | TE internal end (PPT, env) — *conserved across copies* | conserved (LTR family) | random genomic — *non-conserved across copies* |

So if you take all 5'LTRs of one family, extend them by, say, ±50 bp, and align them, you see:

```
                       |---- LTR body ----|
   5' flank             5' boundary       3' boundary             3' flank
   (genomic, random) ---|------ conserved ------|--- TE internal start ---
                                                  (PBS / first 50 bp of internal,
                                                   conserved across all copies)
```

The **5' end of the LTR body** is the position where the per-column conservation jumps from "random" (≈25% per base) to "high" (close to 1.0). The **3' end of the LTR body** is the position where the conservation *stays* high (because the TE-internal start, e.g. PBS, is itself conserved across copies of the same family) — so for 5'LTRs alone, the 3' boundary is **not** marked by a conservation drop.

The 3'LTRs give the mirror picture: the 3' boundary is marked by a conservation drop on the right side, but the 5' boundary is buried in conserved TE-internal-end sequence.

**Combining 5' and 3' LTRs of the same family therefore gives complementary boundary evidence:**

| Boundary to find | Best evidence comes from         | Why                                                   |
|------------------|----------------------------------|-------------------------------------------------------|
| 5' end of LTR    | 5'LTR subset's 5' flank          | Genomic-side flank is non-conserved → sharp transition |
| 3' end of LTR    | 3'LTR subset's 3' flank          | Genomic-side flank is non-conserved → sharp transition |

Inside the LTR body, both subsets are conserved and reinforce each other.

> **Important nuance:** the 3' flank of a 5'LTR is **not** the same sequence as the 5' flank of a 3'LTR. The former is the *start* of the internal region (PBS, gag); the latter is the *end* of the internal region (PPT, env). They are both conserved within their own subset, but they do not align to each other. This means a naive joint MSA of 5'LTRs + 3'LTRs *with flanks* will produce a confused alignment in the flank region. The design below handles this explicitly.




---

## 3. Design proposal

### 3.1 High-level pipeline

```
DANTE_LTR GFF3
  │
  ▼
[A] Element-pair-preserving clustering (per lineage)
  │      • run MMseqs2 on LTR sequences as today
  │      • POST-PROCESS: union the cluster of every 5'LTR
  │        with the cluster of its sibling 3'LTR (same Parent ID)
  │      • result: clusters in which every element is intact
  ▼
[B] Per-cluster flank-aware MSA
  │      • extract each LTR ± F bp flank (default F = 50)
  │      • tag each member with its role: 5LTR or 3LTR
  │      • run MAFFT on the LTR-body region only (anchored MSA),
  │        OR run MAFFT on the full extended sequence and mask flanks
  │        from alignment guidance — see §3.3 for two options
  ▼
[C] Per-subset conservation profiles
  │      • compute per-column conservation separately for the
  │        5LTR and 3LTR subsets
  │      • compute joint conservation across all members
  ▼
[D] Boundary detection
  │      • 5' boundary  ← change-point on 5LTR-subset profile (left of body)
  │      • 3' boundary  ← change-point on 3LTR-subset profile (right of body)
  │      • cross-check that body length is consistent with input annotation
  ▼
[E] Boundary-aware consensus
  │      • build column-wise majority consensus
  │      • emit only the columns between detected 5' and 3' boundaries
  │      • NO 50%-gap filter inside the body (it can erode the ends);
  │        gap handling is replaced by boundary detection
  ▼
[F] (Optional) per-LTR boundary correction
         • for each input LTR, project the cluster's boundaries
           back onto the genome and emit a corrected GFF3 record
         • this is the highest-value downstream improvement: it
           feeds back into solo-LTR detection accuracy
```

### 3.2 Element-pair-preserving clustering (step A)

**Problem.** MMseqs2 may put a 5'LTR and its sibling 3'LTR into different clusters when one of them is degraded or partial. This splits biologically-paired sequences and weakens the per-cluster MSA.

**Proposed algorithm.**

1. Run MMseqs2 within a lineage exactly as today, producing `member → representative` membership.
2. Build a graph `G` on the LTR set:
   - For every MMseqs2 cluster, add edges between all members of that cluster (or at minimum a star to the representative).
   - For every TE element, add an edge between its 5'LTR and its 3'LTR (use `Parent` ID from the GFF3 to identify siblings).
3. Take connected components of `G`. Each component is one final cluster.

This guarantees that 5'/3' siblings are always co-clustered, while still allowing MMseqs2 to merge unrelated copies into family-sized groups. The cost is that two MMseqs2 clusters that share even one element pair become one — which is exactly what we want.

**Edge case:** elements with only one LTR annotated (e.g. the 3'LTR is missing in the input). These contribute a single node, no sibling edge, and rely entirely on MMseqs2 grouping.

### 3.3 Flank-aware MSA construction (step B)

We need an MSA whose **body region** is well-aligned across both 5'LTRs and 3'LTRs, but whose **flank regions** are usable for boundary detection. As noted in §2, naively aligning extended sequences from both subsets is dangerous because the two flank "ends" are different conserved sequences.

**Two viable strategies. Pick one in implementation, design here keeps both open:**

#### Option 1 — Anchored MSA (recommended)

1. Extract each LTR with flanks: `(start - F, end + F)`, F ≈ 50 bp.
2. Build the MSA in two phases:
   - **Phase 1:** align *only* the annotated LTR-body region (no flanks). MAFFT on body sequences only. This is well-conditioned because all members are similar.
   - **Phase 2:** for each member, paste back its real flanks at the appropriate ends of its aligned row, padded with `-` to match alignment columns. The flanks themselves are *not* aligned between members — they sit in unaligned offset columns relative to the body anchor.
3. Compute conservation column-wise in the body-aligned region as usual. For the flank region, conservation must be computed using the *raw* (unaligned) flank columns, indexed by **distance from the body boundary**, not by MSA column.

This is conceptually clean: the body MSA gives positionally-comparable columns; flanks are compared "by offset from the LTR boundary" within each subset.

#### Option 2 — Joint MSA, post-hoc subset conservation

1. Extract each LTR with flanks as above.
2. Run a single MAFFT alignment on all extended sequences. MAFFT will align the LTR body well; the flanks will mostly fall into low-conservation columns (because the two subsets are mutually unrelated there).
3. When computing conservation profiles, **stratify by subset**: compute conservation among 5'LTR members only and among 3'LTR members only, separately. Each subset's conservation is well-defined column-by-column, even where the cross-subset alignment is meaningless.
4. Boundary detection works on the per-subset profiles, ignoring the joint profile in the flank region.

Option 2 is simpler to implement but makes column-indexing more fragile (a column that contains 5'LTR flank residues and 3'LTR body residues is hard to interpret in the joint view). Option 1 is more code but more robust. **Recommendation: Option 1.**

### 3.4 Boundary detection from conservation profiles (step D)

Define a per-column conservation score `c(i)` ∈ [0, 1] — e.g. the frequency of the most common non-gap base, weighted by `(1 - gap_fraction)` so that mostly-gap columns get a low score.

For each subset:
- Build a profile as a function of *distance from the annotated LTR end*. Position 0 = annotated boundary; positions −F … −1 are inside-the-LTR (last F bp of body); positions +1 … +F are flank.

#### 5' boundary (use 5'LTR subset)

- Walk from position −F (deep inside body) outward (towards +F).
- Find the position where conservation transitions from "high" to "background". Use a sliding window (e.g. 5 bp) and require the mean conservation to drop below `T_low ≈ 0.4` while the inner side is above `T_high ≈ 0.7`.
- That position is the corrected 5' boundary.

#### 3' boundary (use 3'LTR subset)

- Mirror image: walk from inside the body towards the 3' flank, find the same kind of high→low transition.

#### Sanity checks

- The corrected boundary should not be more than `±F` from the annotated boundary (we only have flanks of length F).
- The two ends should produce a body length consistent with the LTR family typical length (compare to median annotated length within the cluster; reject corrections that move the length by > 50% — fall back to annotated boundary in that case).
- If the conservation profile shows no clear transition (e.g. the 5' flank is also conserved — unusual but possible if the family inserts into a specific motif), fall back to the annotated boundary and log a warning.

### 3.5 Boundary-aware consensus (step E)

Once 5' and 3' boundaries are known:

- Build the column-wise majority consensus from the body-aligned region (between the two corrected boundaries).
- **Drop the 50%-gap column filter** that exists today (`build_ltr_library.R:78-83`). It was a heuristic for trimming overhangs; with explicit boundary detection it is unnecessary and can erode the ends of well-supported consensuses.
- Inside the body, gaps are still excluded from the majority vote (as today), but the column itself is kept as long as it lies between the boundaries and has at least one non-gap residue.
- Length sanity: if the resulting consensus is shorter than 50 bp (current threshold), fall back to the highest-rank input LTR with annotated coordinates — but log this as a failure of the cluster, not as a normal path.

### 3.6 Per-LTR boundary correction (step F, optional but high-value)

The MSA gives us, for each member LTR, the column offset between its annotated boundary and the cluster's corrected boundary. We can therefore emit a corrected GFF3 record for each input LTR by shifting its `start` / `end` by that offset on the reference. This:

- Improves the input to downstream `detect_solo_ltr.R` BLAST queries (the library sequences are now precisely the LTR, not LTR ± noise).
- Lets the user inspect how often DANTE_LTR's original boundaries needed correcting, which is itself an interesting QC metric.
- Is opt-in (e.g. behind a flag) so it doesn't change current behaviour by default.

This step requires thought about **strand**: corrections must be applied in the orientation of the original feature. Since `getSeq` already reverse-complements minus-strand features before MSA, the offsets are in "canonical orientation" and need to be flipped back when projecting onto a minus-strand feature.

---

## 4. Algorithmic details and parameters

| Parameter             | Default | Notes                                                      |
|-----------------------|---------|------------------------------------------------------------|
| `flank`               | 50 bp   | Currently declared but unused. Activate it.                |
| `min_cluster_size`    | 3       | Keep current default.                                      |
| `conservation_window` | 5 bp    | Sliding window for change-point detection.                 |
| `T_high`              | 0.7     | "Inside body" conservation threshold.                      |
| `T_low`               | 0.4     | "Background" conservation threshold.                       |
| `max_boundary_shift`  | F       | Cannot move boundary further than the available flank.    |
| `max_length_change`   | 50%     | Sanity-check on resulting body length.                     |
| `body_only_first_pass`| true    | Use Option 1 (anchored MSA) by default.                    |

### Conservation score formula (proposed)

For column `i`:
```
gap_frac(i) = #gaps / N
non_gap     = column with gaps removed
top_freq(i) = max freq of any base in non_gap, or 0 if non_gap is empty
c(i)        = top_freq(i) * (1 - gap_frac(i))
```

For the per-subset profile, restrict the column to members of that subset (5LTR or 3LTR) before computing.

### Change-point detection

Simple, robust:
1. Compute `c(i)` for the window of interest.
2. Apply a sliding mean with `conservation_window` width → `cw(i)`.
3. Walk from inside the body outward; find the **first** position where `cw(i) < T_low` *and* there exists a position 1 window earlier where `cw > T_high`.
4. That's the boundary.

If neither transition is found, return "no correction".

Avoid full segmentation libraries — keep the implementation in base R / Biostrings.

---

## 5. Edge cases and open questions

1. **Singleton clusters.** Cannot do boundary detection. Behaviour unchanged from today (use the highest-rank single LTR, no MSA).
2. **Clusters with only 5'LTRs (or only 3'LTRs).** Half the boundary evidence is missing. Detect the available boundary; for the other end, fall back to annotated coordinates.v COMMENT -clustering with 5'LTR only is prefered approach.
3. **Solo LTRs in input.** If the input GFF3 contains pre-existing solo-LTR features, they have no internal flank and would degrade the MSA's flank signal. **Filter them out before clustering** — only paired LTRs (those whose Parent has a complete `transposable_element` with both 5' and 3' annotations) feed into the MSA. - COMMENT - this will not happen and normap pipeline return only pairs of LTRs.
4. **Truncated / partial LTRs.** Acceptable; the MSA majority will out-vote them. The boundary detection uses the *consensus*, not individual members, so a few short members do not move the boundary. COMMENT - aggreed - the majority will outvote the truncated LTRs and the boundary detection is based on the consensus, so it should be robust to a few truncated members.
5. **Highly diverged families.** If MMseqs2 over-splits and the element-pair union doesn't merge enough, conservation in the body may be too low for the boundary thresholds to fire. Tunable, and the fallback (annotated coordinates) is safe.
6. **Insertion-site preference.** Some retrotransposons insert into specific motifs (e.g. tRNA genes). The "genomic" flank of such elements is *not* fully random. This can blur the 5'-boundary detection for the 5'LTR subset. Mitigation: use a stricter `T_low` and require a clear differential, or explicitly check that conservation in the flank is bimodal (mix of conserved + random) versus uniformly conserved. COMMENT - not important at this stage
7. **Element pairs that are actually unrelated.** A mis-annotated TE might pair a 5'LTR with a 3'LTR from a different family. The element-pair union step would then incorrectly merge two clusters. Mitigation: before unioning, check sequence identity between the two siblings (e.g. require ≥ 70% identity) and skip the union if the check fails — log such events as suspect element annotations. COMMENT - gived DANTE_LTR annotation this should not be issue
8. **Reverse-complement consistency.** Already handled by `getSeq` on stranded GRanges, but the new pipeline must keep the role tag (5LTR / 3LTR) attached to the *canonical-orientation* sequence, not the genomic-orientation one. Verify with a test case on a minus-strand element.

---

## 6. Outputs and backwards compatibility

The improved builder produces the same set of files as today:

- `*_LTR_library.fasta` — now boundary-corrected consensuses
- `*_LTR_library_map.tsv` — unchanged structure
- `*_5UTR_tags.fasta`, `*_PPT_tags.fasta` — unchanged (these come from the original annotation, not the consensus)
- `*_tsd_length_map.tsv` — unchanged
- `mafft_alignments/` — now also contains the per-cluster flank-extended MSAs and a small TSV per cluster recording: `(annotated_5_boundary, corrected_5_boundary, annotated_3_boundary, corrected_3_boundary, n_5LTR, n_3LTR, body_length)`

New optional outputs (behind a flag):

- `*_LTR_library_boundary_corrections.tsv` — per-input-LTR offset deltas, useful for QC
- `*_corrected_ltr_features.gff3` — input LTRs with corrected coordinates

No file format changes; the downstream `detect_solo_ltr.R` consumes the same `LTR_library.fasta` and `LTR_library_map.tsv` it does today.

---

## 7. Validation strategy

Before/after comparisons on `test_data/`:

1. **Library sequence length distribution.** Expect tighter distribution after correction; expect mean length close to the lineage-typical LTR length from `databases/lineage_domain_order.csv`.
2. **BLAST self-search.** Each consensus BLASTed against itself should return a single full-length hit; with poor boundaries you sometimes get two staggered hits. Count the fraction of clean self-hits.
3. **Solo-LTR call counts.** Run `dante_ltr_solo` end-to-end with the new library and compare:
   - Number of solo LTRs called
   - Per-lineage `Rsf` (solo / complete ratio) — should be more stable across closely-related lineages
   - Number of solo LTRs whose BLAST hit covers ≥ 95% of the library sequence (should go *up*, since the library is now sharper)
4. **Per-LTR correction sanity.** Distribution of boundary shifts. Most shifts should be ≤ 5 bp; large shifts should be flagged for manual inspection.
5. **Regression on fully-clean elements.** Pick a handful of well-curated elements where boundaries are known to be correct in the input. The pipeline must not move those boundaries.

---

## 8. Out of scope (for this design)

- Re-running DANTE_LTR with corrected boundaries to see if classification ranks change. Possible follow-up.
- Building consensus *across* lineages (e.g. merging closely related Ty3/gypsy sublineages). Today's per-lineage partition is preserved.
- Replacing MAFFT with another aligner. MAFFT `--auto` is adequate for the cluster sizes seen here.
- Changing the 5'UTR / PPT tag database construction. Those are independent of the consensus pipeline.

---

## 9. Summary of the change in one paragraph

Currently the LTR library is a column-wise majority consensus of MMseqs2-clustered LTR features, built from the exact annotated coordinates with no flanks and no boundary refinement. The proposed change adds (a) **element-pair-preserving clustering**, so siblings of the same TE are always grouped together; (b) **flank-aware MSA**, so each LTR is aligned with ± 50 bp of its genomic context; (c) **subset-stratified conservation profiles**, exploiting the asymmetry that 5'LTR genomic flanks are non-conserved on the 5' side while 3'LTR genomic flanks are non-conserved on the 3' side; (d) **change-point boundary detection** on those profiles, which corrects the cluster's LTR boundaries even if every individual input LTR has a slightly wrong annotation; and (e) optionally, **back-projection of corrected boundaries** onto the per-LTR GFF3 records. The result is a sharper library that yields better solo-LTR BLAST hits and exposes annotation errors as a QC by-product.
