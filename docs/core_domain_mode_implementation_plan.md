# Core-domain mode — implementation plan (v1)

**Design:** `docs/core_domain_mode_design.md` (all decisions resolved,
§12; nothing open, §13).
**Scope:** concrete file-level changes against current `main`
(dante_ltr v0.5.4.0).
**Testing discipline:** `smoke` + `short` + the new `core` target must
pass locally before the first commit of each step. Lineage-mode output
must stay **MD5-identical** throughout — that is the gate that says this
feature is additive. No `git push` until the whole plan lands.

---

## 1. Files touched

```
databases/
└── core_domain_order.csv                 [new]  superfamily constraints (design §5)

utils/
├── core_ltr_utils.R                      [new]  merging, seeding, walk, classification
├── detect_core_ltr.R                     [new]  per-chunk driver, mirrors detect_putative_ltr.R
├── calibrate_core_constraints.py         [new]  already written, untracked
└── measure_core_seeding.py               [new]  already written, untracked

dante_ltr                                 [modified]  --mode dispatch + new flags

tests/
├── core.sh                               [new]
├── core_selftest.R                       [new]  unit tests for core_ltr_utils.R
└── data/
    └── core_drapa/                       [new]  73 kb Drapa window, the motivating case

tests.sh                                  [add 'core' dispatcher target]
.github/workflows/tests.yml               [add ./tests.sh core step]
README.md, CLAUDE.md, changelog.md, version.py   [modified]
```

**Explicitly not touched:** `utils/ltr_utils.R`,
`utils/detect_putative_ltr.R`, `dante_ltr_to_library`, `dante_ltr_solo`,
`dante_ltr_summary`, `clean_ltr.R`.

Two facts make that possible and are worth stating up front, because
they are what keeps this plan small:

- `dante_filtering()` (`ltr_utils.R:951`) already takes
  `min_similarity`, `Relative_Length` etc. as **parameters with
  defaults**. The two-tier filter of design §6.0.2 is therefore two
  calls with different arguments — no edit to `ltr_utils.R`.
- `get_te_gff3()` (`ltr_utils.R:724`) needs only
  `list(domain = <data.frame>, ltr_info = <list>)`. Core mode can hand
  it the *post-hoc* in-element domain set (§7) rather than a
  pre-clustered one, because by the time it is called the LTR positions
  are already known.

---

## 2. Implementation order (each step is a separate commit)

### Step 1 — Constraints table and measurement tools

1.1 `databases/core_domain_order.csv`, verbatim from design §5:

```
Superfamily	core_order	offset5prime	offset3prime	core_max_gap	core_max_span	ltr_length	accessory_5	accessory_3
Class_I/LTR/Ty1_copia	INT RT RH	9000	8500	2000	6000	100	PROT GAG	-
Class_I/LTR/Ty3_gypsy	RT RH INT	17000	14000	2000	6000	100	PROT GAG	aRH CHD CHDCR
```

1.2 `git add` the two measurement scripts already written. They are
    referenced by the design (§5.1, §5.3) and one of them is asserted
    against in Step 8.

1.3 Sanity check the table parses with the same idiom
    `detect_putative_ltr.R` uses for `lineage_domain_order.csv`:

```r
t <- read.table("databases/core_domain_order.csv", header = TRUE, sep = "\t")
stopifnot(nrow(t) == 2, all(c("Superfamily", "core_order") %in% names(t)))
```

**Commit 1:** `core: ship superfamily constraints table and calibration tools`

---

### Step 2 — Test fixtures (the local iteration vehicle)

Build these before any R code, so every later step has something fast to
run against.

2.1 Reuse `tests/data/smoke` directly (40 kb slice with one annotated
    DLTP element) rather than copying it to `core_smoke/` — core mode
    must find **the same element with the same boundaries**, and a
    duplicated 40 kb genome in the repo buys nothing. This is the
    concordance gate in miniature and it runs in seconds.

2.2 `tests/data/core_drapa/` — carve a ~200 kb window from
    `/mnt/ceph/454_data/Drapa/hifiasm/assembly_2025_07_30/DRA_2025_07_30/`
    `output/hifiasm_assembly.bp.p_ctg.gfa.fasta`
    plus the matching records from
    `.../output/analysis/repeat_annotation/output/DANTE/DANTE.gff3`.
    Pick a window containing at least two ordered gypsy cores, at least
    one with a CHD immediately 3' of the core (which §5.4 shows is the
    common case and exercises the traversal rule).

    Selection must be scripted and recorded in
    `tests/data/core_drapa/README.md` — source path, contig, coordinates
    — so the fixture is reproducible.

2.3 Assert the premise on the fixture before writing code: lineage mode
    on `core_drapa` produces **zero** non-`D` elements. If it does not,
    the window is unrepresentative; pick another.

**Commit 2:** `core: add test fixtures for core-domain mode`

---

### Step 3 — `core_ltr_utils.R` part 1: merging and seeding

Design §6.0.1, §6.1, §6.2.

3.1 `merge_split_domains(g, max_gap = 50, max_db_overlap = 0.3)`
    — GRanges in, GRanges out. Merge same-`Name`, same-strand adjacent
    features whose `Best_Hit_DB_Pos` intervals tile complementary parts
    of the reference in element order. Merged feature spans both, sums
    `Relat_Length` (capped at 1), takes `max(Similarity)`, and gains
    `Nfragments`.

    Runs **before** `dante_filtering` — this ordering is the whole point
    (design §5.3) and must be asserted in Step 8.

3.2 `core_candidates(g, min_rel_length_core)` — features with
    `Name %in% c("RT","RH","INT")` and LTR support, where LTR support is
    `Final_Classification` starting `Class_I|LTR` **or** any
    `Region_Hits_Classifications` entry doing so. Mark the latter
    `Domain_LTR_Support = "secondary"`.

    Note `Region_Hits_Classifications` is currently discarded at the top
    of the `repeat{}` block in `detect_putative_ltr.R`; core mode must
    read it before that point.

3.3 `find_core_seeds(candidates, constraints)` — per seqname, per
    strand: build runs (consecutive candidates within `core_max_gap`),
    enumerate index triples matching `RT,RH,INT` (gypsy) or `INT,RT,RH`
    (copia) read in element orientation, subject to per-gap and
    `core_max_span` caps, then greedily accept non-overlapping triples.

    **Ordering must be fully deterministic** — sort by
    `(span, -sum(Similarity), seqname, start)`. The final positional key
    is not cosmetic: without it, ties resolve by whatever order
    `expand.grid`/`order()` happens to produce, and the repo has already
    been bitten by order-dependent output
    (`docs/dante_ltr_deterministic_clustering_request.md`). Two runs on
    the same input must be byte-identical.

3.4 Self-test harness `tests/core_selftest.R`, run by `tests/core.sh`.
    It sources `utils/core_ltr_utils.R` and builds synthetic GRanges
    directly, so it needs no genome and runs in under a second —
    which is what makes Steps 3-5 iterable:
    - tandem array `RT RH INT RT RH INT` → two seeds
    - minus strand → same seeds, orientation reversed
    - `RT RH` only → no seed
    - gap exceeding `core_max_gap` → no seed
    - duplicated RT within one run → one seed, leftover RT present
    - split RT with complementary `Best_Hit_DB_Pos` → merged
    - two genuine RT copies with overlapping `Best_Hit_DB_Pos` → not merged
    - split beyond `--split_max_gap` → not merged
    - two sub-threshold fragments → rescued by merge-then-filter order

**Commit 3:** `core: domain fragment merging and ordered core seeding`

---

### Step 4 — `core_ltr_utils.R` part 2: search-space delimitation

Design §6.3. This is the only part with no lineage-mode analogue and
the part §5.4 shows is decisive (99 % of Drapa seeds have a domain
within 500 bp of the core edge).

4.1 `seed_search_limits(seed, g_blocking, constraints)` — walk outward
    from the seed on each side over the **standard-threshold** domain
    set. A domain is transparent only if all six conditions of design
    §6.3 hold; anything else blocks. Return per-seed left/right limits.

    Rule 6 (superfamily agreement) fires **only on positive
    disagreement**; a domain resolving no deeper than `Class_I|LTR`
    passes. Record the count as `Accessory_Superfamily_Mismatch`.

4.2 `core_ranges_left()` / `core_ranges_right()` — reuse the exact
    geometry of `get_ranges_left()` / `get_ranges_right()`
    (`ltr_utils.R:241,254`), including the `+100` over-extension into
    the blocking feature and the `offset2 = 300` overlap into the core,
    substituting the per-seed limits from 4.1 for `upstream_domain` /
    `downstream_domain`. Clamp to `seqlengths`.

4.3 Self-tests (added to the Step 3 harness) — one per row of the §6.3
    table, plus the two negative controls that matter: a GAG resolving
    only to `Class_I|LTR` must **not** block (rule 6 must not fire on
    missing information), and a seed with no accessory domains at all
    must run to the offset cap.

**Commit 4:** `core: blocking-by-default search window delimitation`

---

### Step 5 — `core_ltr_utils.R` part 3: reannotation and classification

Design §7.

5.1 `lca_classification(labels)` — lowest common ancestor over
    `|`-separated REXdb paths.

5.2 `classify_core_element(domains, superfamily, lineage_info,
    max_missing_domains)` returning the attribute set of design §8:
    `Final_Classification` (LCA clipped at the order-derived
    superfamily), `Lineage_Call`, `Lineage_Candidates`,
    `Lineage_Support`, `Classification_Demoted`,
    `Classification_Conflict`.

    Order-compatible lineages come from `order_compatible_with()`, a
    subsequence test against rows of `lineage_domain_order.csv` under
    the called superfamily. It replaces the `domain_distance()` reuse an
    earlier draft assumed: that function recycles when the query carries
    more domains than the reference, which is core mode's normal case.
    See design §7.3.

5.3 Self-tests: the four cases of design §7.6 — unanimous lineage;
    two lineages sharing a parent; cross-superfamily disagreement
    (conflict, order wins); all domains unresolved (report superfamily).

**Commit 5:** `core: post-hoc reannotation and LCA classification`

---

### Step 6 — `utils/detect_core_ltr.R` driver

Mirrors the structure of `detect_putative_ltr.R` — same `optparse`
block, same empty-input guards, same output file set.

6.1 Options: `-g -s -o -c -M -L -d -t`, plus `--min_similarity`,
    `--min_relative_length_core`, `--core_max_gap`, `--core_max_span`,
    `--min_ltr_length`, `--max_te_length`, `--core_require_tsd`,
    `--split_max_gap`.

6.2 Pipeline:

```
import gff3 → CHD_CHDCR_correction → gff_cleanup_overlaps
  → merge_split_domains                                     (§6.0.1)
  → g_block <- dante_filtering(standard thresholds)          (§6.0.2)
  → g_seed  <- dante_filtering(core-relaxed), core names only (§6.0.2)
  → seeds <- find_core_seeds(core_candidates(g_seed))        (§6.2)
  → limits <- seed_search_limits(seeds, g_block)             (§6.3)
  → grL/grR → getSeq → get_TE(...)          [ltr_utils.R, unchanged]
  → structural gates G1-G6                                   (§6.4)
  → in-element domain set from g_block
  → get_te_gff3()                           [ltr_utils.R, unchanged]
  → overwrite Name/Final_Classification from classify_core_element()
  → add_pbs / add_pbs_hemi / get_te_rank    [ltr_utils.R, unchanged]
  → rank-D track via get_domain_clusters_alt + trim_gr        (§6.5)
  → export gff3 + statistics + per-rank fasta
```

   Two details that are easy to get wrong:

   - **`get_TE()` takes the seed, `get_te_gff3()` takes the element.**
     At `get_TE` time the element is not yet delimited, so pass the
     seed's three core domains. Once `ltr_info` comes back the span is
     known; collect every `g_block` domain inside it on the element
     strand and pass *that* set to `get_te_gff3()`. This is why no
     change to `ltr_utils.R` is needed.
   - `get_te_gff3()` sets `TE$Name <- TE$Final_Classification <-
     D$Final_Classification[1]`, i.e. from the first domain. Core mode
     must overwrite both from `classify_core_element()` immediately
     after, on the TE **and** both LTR features.

6.3 `source = "dante_ltr_core"` in column 2. Keep the `TE_%08d` +
    seqname ID scheme identical to lineage mode so the Python
    coordinate-remapping and `get_unique_features()` paths work
    unchanged.

6.4 Statistics CSV must keep its exact shape (rows = classification,
    columns = `D, DL, DLT, DLP, DLTP, RT_domain`) so
    `sum_up_stats_files()` merges chunks unchanged — reuse
    `get_te_statistics()`.

6.5 Emit the per-chunk log block of design §8.

**Commit 6:** `core: per-chunk detector utils/detect_core_ltr.R`

---

### Step 7 — Python `--mode` dispatch

7.1 **First, unify the two command-construction sites.** The chunked
    path builds its command through `_detect_ltr_cmd()` (line ~274) but
    the single-chunk path at line ~1591 duplicates the same list
    inline. Make the single-chunk path call `_detect_ltr_cmd()` too.

    Do this as its own change and verify lineage-mode output is
    MD5-identical *before* adding `--mode`, so the refactor and the
    feature cannot mask each other.

7.2 Add `--mode {lineage,core}` (default `lineage`) and the core-mode
    flags of design §3.1. In `_detect_ltr_cmd()`, select the script
    name from the mode and append core flags only in core mode.

7.3 In core mode `--te_constrains` is read as the *superfamily* table
    (design §5); when unset, default to
    `databases/core_domain_order.csv`.

7.4 `--fallback_mode` needs no special handling — it rewrites
    `gff3_input` before the dispatch — but print a note when both are
    given (design §3.3).

**Commit 7:** `core: --mode flag on dante_ltr`

---

### Step 8 — Tests

8.1 `tests/core.sh`, modelled on `tests/smoke.sh`:
    - run the R self-test harness (Steps 3-5)
    - `--mode core` end-to-end on `core_smoke`; assert the known DLTP
      element is recovered within ±20 bp and its superfamily matches
    - `--mode core` end-to-end on `core_drapa`; assert **≥ 1 non-`D`
      element**, against lineage mode's zero on the same input. This is
      the test that encodes the point of the feature.
    - determinism: run core mode twice, assert MD5-identical GFF3
    - `calibrate_core_constraints.py` on lineage-mode smoke output
      exits cleanly; on a fixture with no DLT/DLTP it exits 1 with the
      "nothing to measure" message
    - **lineage-mode regression:** `--mode lineage` on `tests/data/smoke`
      is MD5-identical to the pre-change output

8.2 `tests.sh`: add the `core` target and include it in `all`.

8.3 `.github/workflows/tests.yml`: add a `./tests.sh core` step next to
    the existing `fallback` step.

**Commit 8:** `core: tests/core.sh and CI wiring`

---

### Step 9 — Documentation and release

9.1 README: a short section per the house style — one-sentence
    motivation, the mechanism (ordered RT/RH/INT core → superfamily),
    the flag table, and a link to the design doc. **It must state the
    ceiling honestly:** core mode is not strictly more sensitive; on
    well-covered genomes 2-25 % of validated elements lack a complete
    filtered core (design §5.2), and its advantage is on genomes REXdb
    does not cover (§5.4).

9.2 CLAUDE.md: add `detect_core_ltr.R` and `core_domain_order.csv` to
    the architecture tables; note the `--mode` dispatch.

9.3 `changelog.md` entry and `version.py` bump to 0.6.0.0 (new user
    facing mode, no breaking change).

**Commit 9:** `docs: document core-domain mode; release 0.6.0.0`

---

## 3. Local verification gates (mandatory)

Per the slow-CI lesson: everything below runs locally before any push.

| gate | command | must show |
|---|---|---|
| lineage mode unchanged | `./tests.sh smoke && ./tests.sh short` | pass, and `ltr.gff3` MD5 equal to pre-change |
| fallback unaffected | `./tests.sh fallback` | pass |
| new mode | `./tests.sh core` | pass |
| determinism | core mode twice on `core_drapa` | identical MD5 |
| full suite | `./tests.sh all` | pass |

The MD5 equality of lineage-mode output is the single most important
gate. Capture the baselines **before** Step 7:

```bash
./dante_ltr -g tests/data/smoke/dante.gff3 -s tests/data/smoke/genome.fasta \
            -o /tmp/base -c 2
md5sum /tmp/base.gff3 > /tmp/lineage_baseline.md5
```

---

## 4. Risks to re-evaluate after implementation

| risk | how it will show | response |
|---|---|---|
| Seeds do not convert to elements | core mode on Drapa yields ≪ 1 984 elements | expected to some degree — measure the conversion rate, and report it in §10 rather than tuning until it looks good |
| `--min_relative_length_core 0.3` adds seeds but no elements | element count flat across the sweep | raise the default to 0.6; the design says explicitly this default is justified only on seed counts so far (§13) |
| Blocking walk too aggressive on domain-dense genomes | many elements with tiny LTR windows and no LTR found | inspect `Accessory_Superfamily_Mismatch` and the truncation log; relax rule 6 first, it is the newest and least-evidenced |
| Large windows make BLAST the bottleneck | core mode much slower than lineage mode per chunk | the offsets are caps, not targets; profile before lowering them, since §5.2 argues generosity is cheap for correctness |
| Core mode worse than lineage on covered genomes | concordance test below 90 % conditional | do not ship as anything but opt-in; investigate before changing the default |

---

## 4b. Latent issues found in `ltr_utils.R` (not fixed here)

Three surfaced while implementing core mode. All are left alone,
because any change to them moves lineage-mode output, and this plan's
primary gate is that it does not.

- **`domain_distance()` (`ltr_utils.R:285`) recycles.** It computes
  `d_query_p == d_reference_p[d_reference_p %in% d_query_p]` without
  checking lengths, so whenever the query carries more domains than the
  reference, R warns and the returned distance is meaningless. Lineage
  mode can reach this whenever a cluster has more domains than its
  lineage's canonical order. Core mode uses `order_compatible_with()`
  instead (design §7.3).
- **`export()` fails on a zero-length GRanges** — "arguments imply
  differing number of rows: 0, 1". `trim_gr()` exports its input, so the
  rank-`D` track dies when no partial cluster has more than one domain.
  Lineage mode would hit this on a chunk where every cluster is a
  singleton. Core mode guards the call site.
- **`get_te_statistics()` (`ltr_utils.R:988`) drops elements** whose
  classification is not among the RT domains' classifications, including
  them in no row and in no `Total`. Harmless in lineage mode, where an
  element's classification is by construction one of those; core mode
  computes the classification, so it needs
  `get_core_te_statistics()`.

Worth fixing in lineage mode separately, with its own before/after
comparison.

---

## 5. Deliberately deferred

Carried over from design §14, restated so they do not creep in:

- Regrouping the rank-`D` track by domain order instead of shared
  classification (§6.5).
- Feeding core-mode elements to `dante_ltr_solo` (O6).
- Any change to `dante_ltr_to_library`, `dante_ltr_summary` or
  `clean_ltr.R`.
- Widening `core_domain_order.csv` with the Darwin runs — that is a
  data change, not a code change, and can land separately.

---

## 6. If something is wrong in practice

The feature is opt-in behind `--mode core`, so the blast radius is
bounded: reverting Step 7 alone disables it while leaving the R code in
place. Steps 3-6 add new files only. There is no migration and no
output-format change for existing users.

---

## 7. Acceptance checklist

- [ ] `--mode lineage` output MD5-identical to v0.5.4.0 on all existing
      test data.
- [ ] `--mode core` on `core_drapa` finds ≥ 1 non-`D` element where
      lineage mode finds 0.
- [ ] Of lineage-mode `DLT`/`DLTP` elements **with a complete filtered
      core**, core mode recovers ≥ 90 % within ±20 bp (design §10).
- [ ] Order-derived superfamily agrees with lineage-mode classification
      on every element found by both (design §5.2 measured 51 875 / 51 875).
- [ ] Core mode is byte-deterministic across repeated runs.
- [ ] Statistics CSV merges across chunks unchanged.
- [ ] `tests/core.sh` passes in CI.
- [ ] README states the core-completeness ceiling explicitly.
