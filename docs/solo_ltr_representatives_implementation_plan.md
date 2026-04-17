# Implementation plan — solo-LTR representatives + overlap dedup + TSD fix

**Goal.** Collapse overlapping solo-LTR hits to one representative per locus, propagate the originating library ID, and widen TSD detection to use a length range rather than a single modal length. Produce a raw GFF3 (all hits) alongside the representative GFF3.

**Based on:** the analysis in the preceding turn (`tmp/analyze_overlaps.R`, `tmp/analyze_tsd.R`, `tmp/analyze_overlaps2.R`).

---

## 1. Decisions locked in

| Aspect | Choice |
|---|---|
| Overlap grouping threshold | reciprocal overlap ≥ **50 %** of the shorter member (ignore strand) |
| Demoted representatives (`boundary_uncertain`, `class_conflict`) | **kept** in the main output, flagged with attributes |
| TSD scan range | **±1 around the per-lineage mode** from the map (clipped to 3..8); fall back to 4..6 when the map has no entry |
| 1-mismatch TSD tolerance | allow for length ≥ **4** (relaxed from current ≥ 5) |
| Two outputs | `solo_ltr_raw.gff3` (all hits, deduped but not collapsed) + `solo_ltr.gff3` (one representative per locus) |
| LibraryID propagation | new `LibraryID` attribute on every solo_LTR (raw + representative) |

---

## 2. Files touched

1. **`utils/solo_ltr_utils.R`** — widen `check_tsd()`, add `LibraryID` plumbing through `make_solo_ltr_gff3()`.
2. **`utils/detect_solo_ltr.R`** — pass `qaccver` from BLAST into `hits_gr` / `make_solo_ltr_gff3()`.
3. **`utils/select_solo_representatives.R`** — new script. Reads raw GFF3, writes representative GFF3.
4. **`dante_ltr_solo`** (Python wrapper) — rename final output to `_raw.gff3`, invoke the new R script to produce `solo_ltr.gff3`.

No changes to `build_ltr_library.R`, no changes to output schema of the *raw* file (except the added `LibraryID`).

---

## 3. TSD fix (`utils/solo_ltr_utils.R` → `check_tsd`)

### Current behaviour (summary)
- `scan_lengths`: single int from map, or `4:6` fallback.
- Iterates longest → shortest, attempts exact; if length ≥ 5 also attempts 1-mismatch.
- Returns on first match; preference is by length, not by exactness.

### New behaviour

```r
# scan_lengths derivation
if (single-int tsd_length was supplied) {
  scan_lengths <- seq(max(3L, tsd_length - 1L), min(8L, tsd_length + 1L))
} else {
  scan_lengths <- 4L:6L
}
```

Two-pass match strategy, so exact always wins over 1-mismatch:

1. **Pass 1 — exact match**, lengths descending.
   - Return the first exact hit (longest exact wins).
2. **Pass 2 — 1-mismatch**, lengths descending, only for `len >= 4`.
   - Return the first hit.
3. Otherwise → `not_found`.

`TSD_length` in the result carries the matched length (now 4, 5, or 6, not fixed at 5).

### Notes
- Range is clipped to `[3, 8]` defensively.
- `tsd_map` contents don't change (still modal per lineage). The widening happens at query time inside `check_tsd()`.

---

## 4. LibraryID plumbing

### `utils/detect_solo_ltr.R`
- `hits_raw$qaccver` is already parsed. Add `LibraryID = hits_raw$qaccver` to the `GRanges()` metadata (line ~195).

### `utils/solo_ltr_utils.R` → `make_solo_ltr_gff3`
- Copy `hit$LibraryID` to the assembled row:
  ```r
  solo$LibraryID <- as.character(hit$LibraryID)
  ```
- Attribute appears right after `Rank` in the output.

Raw GFF3 rows will then carry e.g. `LibraryID=LTR_000123`.

---

## 5. Representative selection — `utils/select_solo_representatives.R` (new)

### CLI

```
select_solo_representatives.R
  -i INPUT_GFF3        # raw solo GFF3 (with LibraryID)
  -o OUTPUT_GFF3       # representative GFF3
  [--overlap_frac 0.5] # reciprocal-overlap threshold (of shorter member)
  [--nest_frac 0.8]    # SL-inside-SL_noTSD threshold for boundary_uncertain
```

### Algorithm

1. **Read GFF3**, split into `solo <- gff[gff$type == "solo_LTR"]` and the rest (TSDs passed through selectively at the end).

2. **Build overlap graph.**
   - `hits <- findOverlaps(solo, solo, ignore.strand=TRUE, drop.self=TRUE, drop.redundant=TRUE)`
   - For each pair `(i, j)`: compute `ov = min(end_i, end_j) - max(start_i, start_j) + 1`; `frac = ov / min(width_i, width_j)`. Keep pairs with `frac >= overlap_frac`.

3. **Connected components (union-find).**
   - Simple in-script union-find over `length(solo)` nodes. For each kept edge, union the two endpoints.
   - Component labels `→` cluster ID.

4. **Per-cluster representative selection.**
   For each cluster `C`:
   - `SL_mask = solo$Rank == "SL"` in C.
   - If any SL: candidate pool = SL; else candidate pool = SL_noTSD.
   - Sort pool by `(-width, -Identity, LibraryID ascending)`; representative = first.
   - `rep_id = solo$ID[representative_index]`.

5. **Demotion flags** for the chosen representative:
   - **`boundary_uncertain`.** True iff the representative is SL and any SL_noTSD member in the cluster is *strictly longer* and contains ≥ `nest_frac` of the representative:
     ```r
     nest_cov = max(
        ov(SL_rep, noTSD_j) / width(SL_rep)
        for each noTSD j in cluster with width[j] > width[rep]
     )
     boundary_uncertain = (nest_cov >= nest_frac)
     ```
   - **`class_conflict`.** True iff any pair in the cluster has `frac >= overlap_frac` *and* different `Final_Classification`.

6. **Aggregate attributes for the representative:**
   - `ClusterSize = length(cluster)`
   - `SupportingHits = comma-joined non-representative LibraryIDs (unique)`; empty for singletons.
   - `boundary_uncertain = "true"` or absent
   - `class_conflict = "true"` or absent

7. **Keep TSD children** (`target_site_duplication` features) whose `Parent` ID equals the representative's ID. Drop all other TSD children.

8. **Write GFF3** preserving order `seqnames → start`, using `rtracklayer::export()`.

### Union-find (inline)

```r
uf_new   <- function(n) seq_len(n)
uf_find  <- function(p, x) { while (p[x] != x) { p[x] <- p[p[x]]; x <- p[x] }; x }
uf_union <- function(p, a, b) {
  ra <- uf_find(p, a); rb <- uf_find(p, b)
  if (ra != rb) p[ra] <- rb
  p
}
```

Good enough for O(N) hits per chunk; no external `igraph` dependency.

### Complexity sanity check
- `findOverlaps(self)` on 3726 hits produces ≈ O(N) pairs (sparse). Union-find is near-linear. Whole pass runs in seconds.

---

## 6. Wiring in the Python wrapper (`dante_ltr_solo`)

### Common change
- Final user-facing outputs become:
  - `solo_ltr.gff3` — representatives only
  - `solo_ltr_raw.gff3` — all hits (post-dedup of exact duplicates)

### `_run_no_chunking`
- Today: copies `chunk_000_solo.gff3` → `final_gff3`.
- New:
  1. Write `solo_ltr_raw.gff3` by copying the chunk output (it's already one-chunk, no merge needed — but still run `dedup_solo_gff3` here too, to drop exact-coord duplicates).
  2. Run `select_solo_representatives.R -i solo_ltr_raw.gff3 -o solo_ltr.gff3`.
  3. Statistics come from the *representative* file (only one SL/SL_noTSD per locus — what the user cares about).

### `_run_chunked`
- Today: `merged_raw.gff3` → `dedup_solo_gff3` → `final_gff3`.
- New:
  1. Dedup produces `solo_ltr_raw.gff3` (renamed target).
  2. `select_solo_representatives.R -i solo_ltr_raw.gff3 -o solo_ltr.gff3`.
  3. Statistics recomputed by R on the representative output (simpler than re-tallying; reuse `get_solo_ltr_statistics()` from `solo_ltr_utils.R`).

### Subprocess call

```python
def select_representatives(tool_path, in_gff, out_gff):
    cmd = [f'{tool_path}/utils/select_solo_representatives.R',
           '-i', in_gff, '-o', out_gff]
    run_subprocess(cmd, "Selecting representative solo LTRs per locus")
```

Statistics file `solo_ltr_statistics.csv` is regenerated from `solo_ltr.gff3`; the per-chunk CSV merge path is dropped (it was counting duplicates anyway). Simpler, and the numbers finally reflect loci instead of hits.

---

## 7. Expected impact on the current test dataset

- Raw output: ~3726 solo_LTRs (same as today).
- Representative output: ~1898 representatives (the union-merged count from the analysis; reciprocal-50 % may drop this a bit lower — maybe 1870 — because one pair with sub-50 % overlap exists).
- `boundary_uncertain` ≈ 20 (the nested-SL cases).
- `class_conflict` = 0 on this dataset; one row in the 4-chunk adjacency case will NOT be flagged because reciprocal overlap is ~0.
- TSD changes: tiny increment in SL count (maybe +3-5 on this genome; 4-bp Tekay, 6-bp SIRE/Tork).
- Every solo_LTR row now carries `LibraryID=LTR_XXXXXX`.

---

## 8. Validation steps after implementation

1. `./utils/select_solo_representatives.R -i tmp/test_solo_new/solo_ltr.gff3 -o tmp/test_sel.gff3` and compare:
   - number of solo_LTRs before vs after
   - count of `boundary_uncertain`
   - count of `class_conflict`
   - spot-check the locus `ctg1492:1416178-1418274` — should yield one SIRE representative at 1416178-1417285 (the SL with TSD=GAAGA), plus one separate Ogre at 1417285-1418274 (the 1 bp touch should not merge).
2. Re-run `dante_ltr_solo` end-to-end; inspect `solo_ltr.gff3` and `solo_ltr_raw.gff3`:
   - every solo_LTR carries `LibraryID=...`
   - TSD lengths include 4 and 6 (expect few but non-zero)
   - `solo_ltr_statistics.csv` counts reflect the representative file, not the raw file
3. Rerun `tests.sh` to confirm nothing in the existing test path regressed.

---

## 9. Out of scope

- Re-running library build with the new boundary-correction (independent change, already done).
- Re-classifying demoted representatives — we flag them, we don't rewrite `Final_Classification`.
- Any change to PBS / junction checks.
