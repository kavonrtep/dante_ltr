# LTR boundary refinement v2 — implementation plan

**Spec / data backing:** `docs/refine_v2_analysis.md`.
**Scope:** rewrite of `dante_ltr_refine` algorithm and GFF3 emission to
implement the *inner-primary policy* with TSD-aware per-side gating.
**Testing discipline:** each step of §3 is a separate commit; smoke +
`./tests.sh refine` must pass locally before each commit. No `git push`
until the whole plan lands.

---

## 1. Files touched

```
utils/
├── refine_boundaries.py        [major changes]
├── parasail_boundary.py        [minor: helpers for inner-pool extraction]
├── tsd_redetect.py             [new: Python port of R TSD detection]
└── refine_mafft_fallback.R     [minor: use inner-pool validation rate]

dante_ltr_refine                [CLI updates]

tests/
├── refine.sh                   [update assertions to v2 schema]
└── data/refine/                [no fixture change; assertions update]

docs/
├── refine_v2_analysis.md       [exists; the spec]
├── refine_v2_implementation_plan.md [this file]
└── refine_v2_gff_schema.md     [new: short attribute-by-attribute reference]

version.py                      [bump → 0.5.0.0 (algorithm change)]
conda/dante_ltr/meta.yaml       [unchanged: reads version.py via regex]
```

No changes to `dante_ltr`, `dante_ltr_to_library`, `dante_ltr_solo`, or
`dante_ltr_summary` are required. Downstream tools that read
`Refinement_Confidence` must be updated only if they switch on the
specific keyword set (see §8 for the rename).

---

## 2. Algorithm changes (semantic summary)

The current pipeline (`utils/refine_boundaries.py:273` `refine_one_cluster`)
runs parasail per side using only same-role anchors and emits the
result whenever it converges. v2 changes:

1. **Two parasail passes per query × side** — once with the outer-edge
   pool (same-role anchors), once with the inner-edge pool
   (opposite-role anchors used at their inner edge). Existing
   `extract_seq_for_side()` already produces the right geometry for
   both pools (verified in the validation experiment).
2. **Inner-primary policy** for picking the boundary coordinate per
   side: take the inner-pool coordinate when `inner_motif_ok = TRUE`;
   otherwise keep DANTE_LTR's original coordinate.
3. **Confidence label per side** based on outer-pool agreement:
   `dual` / `divergent` / `inner_only` / `unrefined`.
4. **TSD redetection** at proposed coordinates, in Python, mirroring
   `evaluate_ltr()` in `utils/ltr_utils.R:651`. Validated bit-for-bit
   against R on `at` and `g2` in §4.4 of the spec.
5. **TSD-loss revert rule (per element).** If the proposed coords on
   *any* side produce a `lost` TSD outcome (original had a TSD, new
   coords find none), revert *all changed sides* on that element to
   DANTE_LTR original. This couples the two sides whenever a TSD
   existed; one-sided changes that destroy the TSD are abandoned.
6. **TSD-gain attribution.** When a refinement produces a TSD where
   none was originally present (`gained` outcome), update the
   `transposable_element` row's `TSD=` attribute and **add two new
   `target_site_duplication` child rows** with the gained coordinates.
   Child rows inherit the parent TE's `Refinement_Confidence`.
7. **MAFFT fallback** still runs for clusters where the **inner pool**
   (not the outer) fails on a high fraction of members, but the
   threshold and exact wiring are deferred (see §7 of the spec).
   For v2, keep current threshold parametrisation; the fallback's
   per-member proposals must also pass the TSD-loss revert rule.
8. **Cluster size threshold** applies to *each role* independently
   (`min_cluster_size` checks N(5'LTR) ≥ N AND N(3'LTR) ≥ N).
9. **Joint-MSA refinement (option C-3)**: out of v2 scope.

---

## 3. Implementation order

Each step below is a separate commit. The order is chosen so each
step leaves the tool in a runnable, test-passing state.

### Step 1 — Add `utils/tsd_redetect.py` (Python TSD port)

A pure library module. No callers yet.

1.1 Port `evaluate_ltr()` TSD logic from `utils/ltr_utils.R:651-707`
    into `utils/tsd_redetect.py`. Single public function:

    ```python
    def detect_tsd(genome, chrom, ltr5_outer, ltr3_outer, strand,
                   max_len=8, min_exact=4, min_fuzzy=5
                   ) -> tuple[int, str, bool]:
        """Returns (length, sequence, fuzzy_flag).
        length=0 indicates not_found."""
    ```

1.2 Bring over the bit-validated implementation in
    `tmp/analyze_pair_and_tsd.py:_detect_tsd_biological` (cleaned up
    — single function, no dead-code branch).

1.3 Add a tiny self-test at the bottom (run with
    `python utils/tsd_redetect.py --selftest`) that exercises one +
    strand case and one - strand case.

**Commit 1:** `refine: tsd_redetect — Python port of R TSD detection`.

### Step 2 — Refactor `parasail_boundary.py` extraction helpers

`extract_seq_for_side()` already supports both pools (verified in
validation script). Two small additions:

2.1 Add convenience function:

    ```python
    def extract_pool_for_side(members, role, side, genome, anchor_len,
                               flank_len, genome_lens) -> list[dict]:
        """Filter `members` to `role`, extract per-member seq for `side`.
        Returns the list expected by aggregate_extension."""
    ```

2.2 Extract `aggregate_extension()` from `process_cluster_side` into a
    standalone function that takes a query-info dict and a list of
    anchor-info dicts (matches the validation script's design — same
    signature). Keep `process_cluster_side()` as a thin wrapper
    calling `aggregate_extension()` with same-role anchors, for
    backward compat during the refactor.

**Commit 2:** `refine: parasail_boundary — factor out aggregate_extension and pool helpers`.

### Step 3 — Add `RefRecord` fields and the per-side dual-pool computation

In `utils/refine_boundaries.py`:

3.1 Extend `RefRecord` with new fields:

    ```python
    # Outer-pool result (legacy "A")
    outer_corrected_g: Optional[int] = None
    outer_motif_ok: Optional[bool] = None
    outer_n_pairs: int = 0
    # Inner-pool result (legacy "B")
    inner_corrected_g: Optional[int] = None
    inner_motif_ok: Optional[bool] = None
    inner_n_pairs: int = 0
    # Original-state capture (for the TSD/motif gates)
    orig_motif_ok: Optional[bool] = None
    orig_tsd_len: int = 0
    orig_tsd_seq: str = ""
    # Final after gates (renamed from final_*; semantics changed)
    final_corrected_g: Optional[int] = None
    confidence: str = "unrefined"           # dual / divergent / inner_only / unrefined
    refinement_method: str = "none"         # parasail_inner / parasail_outer / mafft / none
    revert_reason: str = ""                 # "" / "tsd_lost" / "motif_fail" / ...
    tsd_outcome: str = ""                   # kept_exact / shifted / lost / gained / both_none
    tsd_new_len: int = 0
    tsd_new_seq: str = ""
    ```

    Drop fields that no longer apply: `parasail_motif_at`,
    `parasail_snap_offset` (move snap_offset under each pool),
    `final_motif_ok` (replaced by `confidence`),
    `mafft_motif_ok` / `mafft_corrected_g` (move under explicit
    method).

3.2 Replace `refine_one_cluster()` with `refine_one_cluster_v2()`:

    ```python
    def refine_one_cluster_v2(cluster_members, cluster_id, lineage,
                              genome, genome_lens, opts):
        # 0. Capture original motif and TSD state (per side / per element)
        # 1. For each side ('5', '3'):
        #      build outer pool (same role) and inner pool (opposite role)
        #      for each query (same-role member of side):
        #          aggregate_extension(query, outer_anchors)  → outer_g
        #          aggregate_extension(query, inner_anchors)  → inner_g
        #          snap each to motif (TG/CA, ±5 bp)
        # 2. Apply inner-primary policy per side → proposed_g
        # 3. Assign confidence label per side
        # 4. Per-element TSD recheck at proposed coords
        # 5. If TSD outcome == "lost" AND original had TSD: revert all
        #    changed sides on this element
        # 6. Build RefRecord per LTR row, attach TSD outcome + new TSD
        return records
    ```

3.3 Drop the cluster-size filter to require **both** roles ≥ N (was
    `len(cluster) ≥ N`). Document the change in `--min_cluster_size`
    help text.

3.4 Keep `cluster_validation_rate()` but compute it against
    `inner_motif_ok` (was `parasail_motif_ok`). This is what feeds
    the MAFFT fallback decision.

**Commit 3:** `refine: dual-pool parasail with inner-primary policy + TSD gate`.

After this commit the tool runs end-to-end; the GFF emission still
uses old attributes (Step 4).

### Step 4 — Update GFF3 emission (`emit_refined_gff3`)

The two big changes: TSD-gain inserts new child rows; attribute names
change.

4.1 Output GFF3 attribute schema (specified in
    `docs/refine_v2_gff_schema.md`):

    On `transposable_element` and `long_terminal_repeat` rows:

    | attribute | values |
    |---|---|
    | `Refinement_Method` | `parasail_inner` / `parasail_outer` (rare) / `mafft` / `none` |
    | `Refinement_Confidence` | `dual` / `divergent` / `inner_only` / `unrefined` |
    | `Cluster_ID`, `Cluster_Size` | as before |
    | `Original_Start`, `Original_End` | only when boundary changed |
    | `Motif_Orig`, `Motif_New` | `TG` / `CA` / `none` (the 2 bp at each coord) |
    | `Outer_Pool_g`, `Outer_Pool_Motif_OK` | per-side outer-pool diagnostics (TE-row carries 5LTR/3LTR pair, comma-separated; LTR-row carries scalar) |
    | `Inner_Pool_g`, `Inner_Pool_Motif_OK` | as above |
    | `TSD_Outcome` | `kept_exact` / `kept_fuzzy` / `shifted` / `lost` / `gained` / `both_none` (TE-row only) |
    | `Revert_Reason` | non-empty only if the side reverted |

    On `target_site_duplication` rows (if present):

    | attribute | values |
    |---|---|
    | `Parent` | TE ID |
    | `TSD_Status` | `original` / `re_identified` / `gained` |
    | `Refinement_Confidence` | inherited from parent TE |

    Drop `TG_OK` / `CA_OK` (legacy; replaced by the explicit
    `Motif_*` and per-pool diagnostics).

4.2 In `emit_refined_gff3()`:
    - Change parsing to also collect input `target_site_duplication`
      rows by parent TE id.
    - For each TE:
      - If `tsd_outcome == "lost"` and revert is in effect → emit
        original TSD rows verbatim.
      - If `tsd_outcome == "kept_exact"` / `kept_fuzzy` / `shifted`
        → emit fresh `target_site_duplication` rows at the new
        flanks (computed from `tsd_new_*`); drop the input ones.
      - If `tsd_outcome == "gained"` → emit two new
        `target_site_duplication` rows; remove any input ones (the
        input had `not_found`, so probably none existed, but be
        defensive).
      - If `tsd_outcome == "both_none"` → drop any input rows; emit
        none.
    - Update the TE-row `TSD=` attribute accordingly: replace
      `TSD=not_found` with the new sequence on `gained`; replace
      with new sequence on `shifted`; keep on `kept_*`; replace with
      `not_found` on `lost`.

4.3 Sanity script: a one-shot Python check (in `tmp/`, not committed)
    that reads the new `_refined.gff3`, reconstructs the
    `transposable_element` ↔ `target_site_duplication` parent links,
    and verifies coordinate consistency (TSD's start within
    `[5'LTR.start - 8 .. 5'LTR.start - 1]` etc).

**Commit 4:** `refine: GFF3 schema v2 — TSD child-row updates + new attribute set`.

### Step 5 — Update `_per_element.tsv` and `_clusters.tsv` schemas

5.1 `_per_element.tsv` columns (add / rename):

    | column | source |
    |---|---|
    | `chrom`, `start_orig`, `end_orig`, `strand`, `role`, `parent_te_id`, `lineage_full` | as before |
    | `cluster_id`, `cluster_size` | as before |
    | `outer_corrected_g`, `outer_motif_ok`, `outer_n_pairs` | new per-pool diagnostics |
    | `inner_corrected_g`, `inner_motif_ok`, `inner_n_pairs` | new |
    | `final_corrected_g` | unchanged in spirit |
    | `refinement_method` | new values: `parasail_inner` / `parasail_outer` / `mafft` / `none` |
    | `confidence` | replaces `final_confidence` with new label set |
    | `revert_reason` | new |
    | `orig_tsd_len`, `orig_tsd_seq` | new |
    | `tsd_outcome`, `tsd_new_len`, `tsd_new_seq` | new |
    | `shift_bp` | unchanged |

5.2 `_clusters.tsv`: keep columns, but the `validation_rate` column
    now corresponds to `inner_motif_ok` (rename to
    `inner_validation_rate`). Add `outer_validation_rate` for
    diagnostics.

5.3 `_run.json`: add a top-level `policy: "inner_primary"` and
    confidence-label histogram (counts of dual / divergent /
    inner_only / unrefined).

**Commit 5:** `refine: TSV/JSON schema v2 — per-pool columns + confidence labels`.

### Step 6 — Update CLI (`dante_ltr_refine`)

6.1 In the argparse setup at `utils/refine_boundaries.py:938`:
    - **Drop** `--boundary_motif` `choices` to just `["TG/CA"]` (or
      remove the option entirely; TG/CA always on). Spec calls for
      forced TG/CA; `none` is no longer offered.
    - Add `--no_tsd_revert` (default off) — disables the TSD-loss
      revert rule. Useful for diagnostics.
    - Document `--min_cluster_size` change (per-role, both sides ≥ N).
    - Help text updated to mention inner-primary policy and the new
      confidence label set.

6.2 Update `dante_ltr_refine --help` smoke-test asserts in
    `conda/dante_ltr/meta.yaml:50` — already a `grep -q usage`; no
    change needed.

**Commit 6:** `refine: CLI — drop motif=none, add --no_tsd_revert`.

### Step 7 — Update `tests/refine.sh`

7.1 Replace assertions:
    - `Refinement_Method=parasail` → `Refinement_Method=parasail_inner`
      (or any of the new keywords)
    - `final_confidence` keyword set: `high|medium|low` →
      `dual|divergent|inner_only|unrefined`
    - `TG_OK=` grep → `Motif_New=` grep
    - Add: count of `TSD_Outcome=lost` rows must be **0** in the v2
      output for the fixture (because the revert rule eliminates
      them).
    - Add: count of `target_site_duplication` rows in output ≥ count
      in input (because of `gained` cases).

7.2 Run on `tests/data/refine` and update any column-index references
    in `awk` (the per-element TSV columns shifted).

7.3 The hybrid (parasail + MAFFT fallback) sub-stage stays; assertions
    update similarly.

**Commit 7:** `refine: tests — refresh assertions for v2 schema`.

### Step 8 — Update `utils/refine_mafft_fallback.R`

Minimal changes to keep the fallback aligned with the new gating:

8.1 The R script reads the per-cluster FASTA (LTR ± mafft_flank) and
    proposes per-member coordinates. v2 uses the same I/O contract.
    Add a new column `mafft_motif_ok` (already present) — if the
    Python caller will run TSD recheck on MAFFT's proposals, no R
    change is needed there.

8.2 Python side: after MAFFT proposes a coordinate, run the same
    motif + TSD-loss gate before accepting. This is wired in
    `run_mafft_fallback()` in `utils/refine_boundaries.py:355`.

**Commit 8:** `refine: MAFFT fallback — gate proposals on motif + TSD-loss`.

### Step 9 — Documentation + version bump

9.1 Write `docs/refine_v2_gff_schema.md` (a 1-page table of all v2
    GFF attributes; pure reference for downstream consumers).

9.2 Update `README.md` boundary-refinement section briefly per the
    feedback memory (one-sentence motivation + mechanism + table +
    link to spec).

9.3 Bump `version.py` to `0.5.0.0` (algorithm change).

**Commit 9:** `release 0.5.0.0 — refine v2 (inner-primary policy)`.

---

## 4. CLI changes — summary

| flag | v1 | v2 |
|---|---|---|
| `--boundary_motif` | `TG/CA` (default) / `none` | dropped; TG/CA always on |
| `--min_cluster_size` | total cluster size ≥ N | N(5'LTR) ≥ N AND N(3'LTR) ≥ N |
| `--no_tsd_revert` | — | new: disables TSD-loss revert rule |
| `--no-mafft-fallback` | unchanged |  |
| all numeric tuning flags | unchanged |  |

---

## 5. GFF3 attribute schema — summary

Detailed in `docs/refine_v2_gff_schema.md`. Quick view:

```
##gff-version 3

chr1 dante  transposable_element 1000 6500 . + . \
    ID=TE_1;Final_Classification=…;TSD=GATCG;Refinement_Method=parasail_inner;\
    Refinement_Confidence=dual;Cluster_ID=…;Cluster_Size=12;\
    Original_Start=1003;Original_End=6498;Motif_Orig=TG/CA;Motif_New=TG/CA;\
    Outer_Pool_g=998,6502;Outer_Pool_Motif_OK=TRUE,FALSE;\
    Inner_Pool_g=1000,6500;Inner_Pool_Motif_OK=TRUE,TRUE;\
    TSD_Outcome=kept_exact

chr1 dante  long_terminal_repeat 1000 1500 . + . \
    Parent=TE_1;LTR=5LTR;Refinement_Method=parasail_inner;\
    Refinement_Confidence=dual;Original_Start=1003;Motif_Orig=TG;Motif_New=TG;\
    Outer_Pool_g=998;Outer_Pool_Motif_OK=TRUE;Inner_Pool_g=1000;Inner_Pool_Motif_OK=TRUE

chr1 dante  target_site_duplication 995 999 . + . \
    Parent=TE_1;TSD_Status=re_identified;Refinement_Confidence=dual

chr1 dante  target_site_duplication 6501 6505 . + . \
    Parent=TE_1;TSD_Status=re_identified;Refinement_Confidence=dual
```

For TE-row attributes that have a `5LTR,3LTR` flavour
(`Outer_Pool_g`, `Outer_Pool_Motif_OK`, `Inner_Pool_g`,
`Inner_Pool_Motif_OK`), the order is always 5'LTR first.

---

## 6. Per-element TSV schema — full column list

```
chrom              str
start_orig         int
end_orig           int
strand             '+' | '-'
role               '5LTR' | '3LTR'
parent_te_id       str
lineage_full       str
cluster_id         str
cluster_size       int
outer_corrected_g  int|NA
outer_motif_ok     TRUE|FALSE|NA
outer_n_pairs      int
inner_corrected_g  int|NA
inner_motif_ok     TRUE|FALSE|NA
inner_n_pairs      int
final_corrected_g  int|NA
refinement_method  parasail_inner|parasail_outer|mafft|none
confidence         dual|divergent|inner_only|unrefined
revert_reason      ""|tsd_lost|motif_fail|...
orig_tsd_len       int
orig_tsd_seq       str
tsd_outcome        kept_exact|kept_fuzzy|shifted|lost|gained|both_none
tsd_new_len        int
tsd_new_seq        str
shift_bp           int|NA
```

---

## 7. Testing

### 7.1 Local smoke (after each commit)

```bash
./tests.sh refine
```

Must pass before each commit. The fixture is small (~1 Mb, 78 TEs, 4
lineages) so this completes in under a minute.

### 7.2 Full validation on `at` and `g2` (after Commit 4 or later)

Re-run the production command on both datasets:

```bash
./dante_ltr_refine -g test_data/at/Alyr_dante_ltr.gff3 \
                   -s test_data/at/Alyr.fasta \
                   -o tmp/refine_v2_at \
                   --threads 6
./dante_ltr_refine -g test_data/g2/DANTE_LTR.gff3 \
                   -s test_data/g2/genome.fasta \
                   -o tmp/refine_v2_g2 \
                   --threads 6
```

Acceptance criteria:

- **TG/CA preservation**: `awk` over `_refined.gff3` for
  `Motif_New=TG/CA` should approach 100 % on `transposable_element`
  rows (allowing for a tiny minority where the original had non-TG/CA
  termini, which DANTE_LTR doesn't normally emit but defensive code
  must handle).
- **TSD outcome distribution**: `lost` should be ≤ 0.5 % (revert rule
  in effect); `kept_exact` should dominate.
- **Confidence label histogram** matches the §4.2 simulation
  (within ±2 percentage points): on `g2`, ~64 % `dual`, ~16 %
  `divergent`, ~10 % `inner_only`, ~10 % `unrefined`.
- **Pair symmetry**: per-element width-asymmetry change ≤ 5 bp on
  ≥ 95 % of elements (the TSD-revert rule should suppress most of the
  asymmetric one-sided changes from the v1-policy simulation).

### 7.3 Regression guard

Run `dante_ltr_to_library --refined_gff3` and `dante_ltr_solo
--refined_gff3` against the v2 output to ensure downstream tools
still parse it (Cluster_ID, Refinement_Method,
Refinement_Confidence are read by these tools — must check).

---

## 8. Migration / compatibility notes

- The output GFF3 is a strict superset of v1's structure (all rows
  still present; `target_site_duplication` rows may be added or
  removed; attributes change).
- Downstream code that reads `Refinement_Confidence == "high"|"medium"`
  must be updated to map to `dual|divergent|inner_only`. Search:
  `grep -rn "Refinement_Confidence" utils/ *.R *.py docs/` to find
  all readers; these likely include `utils/build_ltr_library.R`,
  `dante_ltr_solo`, `dante_ltr_summary`. Audit each before Commit 4.
- `TG_OK` / `CA_OK` attributes are removed. If any downstream tool
  reads them, it must be moved to `Motif_New=` parsing or
  `Inner_Pool_Motif_OK=` parsing.
- The validation TSV column rename (`A_*` / `B_*` / `C_*` are
  *internal-only*; the production TSV uses `outer_*` / `inner_*`
  names per §6).

---

## 9. Risks

1. **MAFFT fallback wiring is deferred.** The fallback may rarely
   contribute under v2 (inner pool's high pass rate makes most
   clusters skip it). If after release we see clusters where the
   inner pool also fails systematically, re-tune the trigger
   threshold. Tracked as the open follow-up in the spec.
2. **Pair symmetry for elements with originally-no-TSD.** The TSD-loss
   rule does not apply to `both_none` elements, so one-sided changes
   on those can still grow asymmetry. The validation showed this is
   uncommon (~5–6 % of elements) and we accept it for v2; revisit if
   downstream tools complain.
3. **Schema rename impact.** Several downstream tools probably read
   the old `Refinement_Confidence` keywords. Step 5 (TSV) and Step 4
   (GFF) bracket the breakage; the audit in §8 must precede Commit 4.
4. **Performance.** Two parasail passes per query (outer + inner pool)
   roughly doubles the parasail compute. On `at` the validation
   script took ~8 minutes vs. ~4 for v1 production. Acceptable but
   noticeable on large genomes; revisit if it becomes a bottleneck.

---

## 10. Open follow-ups carried over from the spec

- MAFFT fallback threshold and per-proposal TSD gate detail.

(Joint-MSA refinement, confidence-label keywords, TSD-gain attribution
are all decided in §6 of the spec and incorporated above.)
