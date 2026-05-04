# Parasail-based hybrid LTR boundary refinement — implementation plan

**Status:** proposed
**Author:** Petr Novak (via agent)
**Date:** 2026-05-04
**Related:**
- `docs/boundary_validation_implementation_plan.md` (TG/CA / TSD per-cluster QC; landed)
- `tmp/parasail_compare/*` (exploratory prototypes + benchmarks; not committed)
- `https://github.com/kavonrtep/CARP` `scripts/global_local_aln.py` (reference algorithm)

## 1. Goal

Replace the current MAFFT-based per-cluster boundary correction in
`build_ltr_library.R` with a **per-element boundary refinement** pipeline
based on **parasail anchored extension**, with **MAFFT change-point as
selective fallback**. Produce two outputs:

1. **Refined DANTE_LTR GFF3** — same schema as input but with corrected
   per-element boundaries and refinement-confidence attributes.
2. **Solo-LTR GFF3** (existing) — produced from a refined library built
   only from validated members.

The hybrid is dramatically faster than full MAFFT refinement (5–30×
depending on genome) and significantly more accurate at the per-member
TG/CA level (Alyr: ~76 % member TG/CA-OK with parasail vs ~3–24 % with
MAFFT change-point at the per-cluster boundary).

## 2. Non-goals

- Not changing DANTE_LTR detection itself; refinement consumes
  DANTE_LTR's output.
- Not introducing a new MSA tool. Parasail is pairwise; MAFFT remains
  the fallback.
- Not removing the existing MAFFT-based boundary correction in v1; it
  becomes the fallback path. Removal can be a later cleanup once
  parasail is stable.
- Not redefining "solo LTR" — solo-LTR detection algorithm is unchanged
  except that it consumes a cleaner library.

## 3. Decisions taken (from design discussion)

The following questions were answered during design:

| # | decision |
|---|---|
| Library composition | **C** — per-cluster representative consensus built from **validated members only**. Coverage retained on under-validated clusters by including a fallback consensus flagged `low_confidence`. Threshold configurable via `--min_validated_members` (default 4). |
| Validation criterion | **TG/CA strict by default** (both motifs must match at corrected per-member boundary). **Configurable** via a global `--boundary_motif` flag and via a future per-lineage column in `lineage_domain_order.csv`. |
| Hybrid policy | **Y** — parasail-first; MAFFT change-point fallback only when parasail's per-cluster validation rate is below `--mafft_fallback_threshold` (default 0.5). |
| Pipeline integration | **Separate `dante_ltr_refine` command**; refined GFF3 is its output. Downstream `dante_ltr_solo` consumes the refined GFF3 (does NOT silently fold refinement in). |
| Failed-cluster handling | **Include in output** with `refinement_confidence=unrefined`; do not drop. |
| Backward compatibility | The current MAFFT-based per-cluster correction in `build_ltr_library.R` **stays alive in v1** (used when `--refined_gff3` is not supplied). Removal deferred to a later cleanup release once parasail-based refinement is the established default. |

## 4. Configuration model for boundary motif

The motif (TG/CA) is the validation target. Some lineages don't follow
it. Configuration is layered:

### 4.1 Global default

```
--boundary_motif TG/CA      # default, validates as 5'=TG, 3'=CA
--boundary_motif none       # validate purely on parasail
                              alignment evidence (cumulative score)
--boundary_motif TG-AC      # rare; placeholder for non-TG/CA families
                              (future)
```

`TG/CA` and `none` are the v1 supported values. Other motifs deferred to
when needed.

### 4.2 Per-lineage override (planned, not v1)

Add a new column `boundary_motif` to
`databases/lineage_domain_order.csv`:

```
Lineage                                 Domains   ...   boundary_motif
Class_I/LTR/Ty1_copia/Ale               GAG ...   ...   TG/CA
Class_I/LTR/Ty3_gypsy/non-chromo/Ogre   GAG ...   ...   TG/CA
Class_I/LTR/<some-future-lineage>       GAG ...   ...   none
```

When the column is missing or empty, fall back to the global value. v1
silently uses the global default and ignores the column; v1.x lights it
up.

### 4.3 Implementation note

Boundary-motif evaluation is centralised in a single function
`validate_motif(seq, side, motif)` so that switching from `TG/CA` to
`none` (or to alternatives) is a one-line policy change.

## 5. Architecture

### 5.1 Pipeline

```
DANTE_LTR.gff3  +  genome.fasta
        │
        ▼
┌──────────────────────────────┐
│  dante_ltr_refine            │
│  utils/refine_boundaries.py  │
│  ├─ parse GFF3, getSeq       │
│  ├─ MMseqs2 cluster (lin)    │
│  ├─ per-cluster:             │
│  │  1. parasail anchored ext │
│  │  2. per-member TG/CA snap │
│  │  3. validation policy     │
│  │  4. MAFFT fallback if low │
│  ├─ emit refined GFF3        │
│  ├─ emit per-element TSV     │
│  └─ emit cluster manifest    │
└──────────────────────────────┘
        │
        ├─ refined.gff3 (Output A)
        ├─ refined_per_element.tsv
        └─ clusters_manifest.tsv
        │
        ▼
┌──────────────────────────────┐
│  build_ltr_library.R         │   (extended)
│  ├─ accepts --refined_gff3   │
│  ├─ build consensus from     │
│  │  VALIDATED members only   │
│  └─ flag low-confidence      │
│     clusters                 │
└──────────────────────────────┘
        │
        ├─ <out>_LTR_library.fasta
        ├─ <out>_LTR_library_map.tsv
        ├─ <out>_5UTR_tags.fasta
        └─ <out>_PPT_tags.fasta
        │
        ▼
┌──────────────────────────────┐
│  dante_ltr_solo              │
│  (existing solo-LTR pipeline)│
│  → uses refined library      │
└──────────────────────────────┘
        │
        └─ solo_ltr.gff3 (Output B)
```

### 5.2 Modules

| module | role |
|---|---|
| `utils/parasail_boundary.py` | core parasail anchored-extension + TG/CA snap. Library-style; importable from refine_boundaries.py and any later tool. Ports CARP's `global_local_aln.py` helpers (per-column scoring, optimal-alignment extraction). |
| `utils/refine_boundaries.py` | top-level refinement engine. Reads GFF3 + genome, runs clustering, runs parasail (and MAFFT fallback), emits refined GFF3 + per-element TSV. |
| `utils/refine_mafft_fallback.R` | MAFFT change-point fallback callable per cluster from refine_boundaries.py via `Rscript` subprocess. Reuses logic in `utils/build_ltr_library.R`. |
| `utils/build_ltr_library.R` | extended to accept `--refined_gff3` and build consensus from validated members. |
| `dante_ltr_refine` | Python CLI wrapper that drives `refine_boundaries.py`, mirroring the style of the existing `dante_ltr` and `dante_ltr_solo` scripts. |
| `dante_ltr_solo` | extended to accept `--refined_gff3`; falls back to `--gff3` for backward compat. |

### 5.3 Data flow per cluster (the core loop)

For each MMseqs2 cluster (per lineage, n_members ≥ min_cluster_size):

```
A. parasail anchored extension
   - extract per-member sequences:
       seq_5 = upstream_flank(F) + first ANCHOR_LEN bp of LTR
       seq_3 = last  ANCHOR_LEN bp of LTR + downstream_flank(F)
   - for each pair (i,j), parasail semi-global with one end fixed:
       side=5: --end 3 (anchor at 3' fixed; 5' free)
       side=3: --end 5 (anchor at 5' fixed; 3' free)
   - aggregate per-member: take Nth-largest extension length
   - snap corrected boundary to nearest TG (5') / CA (3')
     within ±SNAP_WINDOW bp on the member's own raw sequence

B. validation (per element, per side)
   - boundary motif matches per --boundary_motif policy
       (TG at 5', CA at 3') -> validated_5, validated_3
   - or, with --boundary_motif=none, accept any parasail-resolved
     boundary

C. cluster-level validation rate
   = fraction of cluster members where the corrected boundary is
     validated on its respective side (5' for 5'LTR, 3' for 3'LTR)

D. MAFFT fallback (only if rate < --mafft_fallback_threshold)
   - run mafft + change-point on the cluster (via existing
     build_ltr_library.R logic, called as subprocess)
   - for each element where parasail failed validation, replace the
     boundary with MAFFT's per-member position (computed from MSA col +
     cum_mat) IF MAFFT's call validates better
   - record refinement_method=mafft for those rows

E. emit per-element refinement record
```

## 6. CLI

### 6.1 `dante_ltr_refine` (new)

```
dante_ltr_refine
  -g, --gff3 GFF3            DANTE_LTR GFF3 input
  -s, --genome FASTA         Genome FASTA
  -o, --output PREFIX        Output prefix
  --identity 0.9             MMseqs2 cluster identity
  --min_cluster_size 6
  --anchor_len 50
  --flank_len 1000           Cap; per-cluster proportional rule applies
  --snap_window 5            ±bp around parasail boundary to snap to motif
  --boundary_motif TG/CA     {TG/CA,none}  default TG/CA
  --mafft_fallback_threshold 0.5
                              cluster-level validation rate trigger
  --no-mafft-fallback        disable fallback; parasail-only
  --threads 4
  --workers 4                parallel cluster workers
  -v, --verbose
```

Outputs at `<output>_*`:
- `<output>_refined.gff3`       — Output A
- `<output>_per_element.tsv`    — per-element refinement record
- `<output>_clusters.tsv`       — cluster manifest with validation stats
- `<output>_run.json`           — parameters + timing summary

### 6.2 `dante_ltr_solo` (extended)

Add input:
```
--refined_gff3 PATH    use refined GFF3 (preferred)
--gff3 PATH            legacy path; falls back to original DANTE_LTR
                       boundaries when --refined_gff3 not given
```

Solo-LTR output schema is unchanged.

### 6.3 `build_ltr_library.R` (extended)

Add input:
```
--refined_gff3 PATH    if given, build consensus from VALIDATED
                       members only (per --boundary_motif); flag
                       clusters with < min_validated_members as
                       low_confidence
--min_validated_members 4
                       at least this many validated members required
                       to build a "high" cluster consensus; below
                       this, build best-effort consensus and flag
                       low_confidence
```

## 7. Output schemas

### 7.1 Refined GFF3 (Output A)

Same column structure as DANTE_LTR.gff3.

For features that were **refined**, modified attributes:

```
transposable_element  start=<refined>  end=<refined>
  attributes:
    ID=<unchanged>
    Final_Classification=<unchanged>
    Rank=<unchanged>
    flanking_5end=<unchanged>
    flanking_3end=<unchanged>
    Original_Start=<original_start>
    Original_End=<original_end>
    Refinement_Method=<parasail|mafft|none>
    Refinement_Confidence=<high|medium|low|unrefined>
    Cluster_ID=<MMseqs2_cluster_rep>
    Cluster_Size=<n>
    TG_OK=<TRUE|FALSE|NA>      (only meaningful when motif=TG/CA)
    CA_OK=<TRUE|FALSE|NA>
```

For unrefined features (parasail and MAFFT both unable to call), keep
original `start`/`end` and add only `Refinement_Method=none` /
`Refinement_Confidence=unrefined` so consumers can filter.

`long_terminal_repeat` features get the same treatment.

### 7.2 Per-element TSV

```
ltr_id   parent_te_id   chrom   start_orig   end_orig
strand   role           lineage_full
cluster_rep            cluster_size
parasail_n_pairs       parasail_corrected_pos      parasail_tg_ok  parasail_ca_ok
mafft_corrected_pos    mafft_tg_ok                 mafft_ca_ok
final_corrected_pos    final_method                final_confidence
shift_bp                                            (final_corrected_pos - original)
```

One row per (LTR, side). 5'LTRs have side=5 (TG), 3'LTRs have side=3
(CA); when motif is TG/CA. With motif=none the boolean columns are NA.

### 7.3 Cluster manifest

```
cluster_id   lineage_full   n_total   n_5ltr   n_3ltr
n_validated_5  n_validated_3
parasail_validation_rate    mafft_invoked    mafft_validation_rate
final_consensus_built_from  (validated|all_members)
low_confidence              (TRUE|FALSE)
```

### 7.4 Library FASTA + map

Same format as today; cluster_id stays compatible. Map TSV gains:

```
ltr_id   Final_Classification   n_validated   n_total   consensus_built_from
                                                         low_confidence
```

## 8. Implementation milestones

Numbering for tracking; each lands as one or a small number of commits.

### M1 — Library-ize parasail (3–4 days)

- Extract working code from `tmp/parasail_compare/run_parasail_boundary.py`
  into `utils/parasail_boundary.py` as importable functions:
  - `parsail_align_pair(seq1, seq2, end, ...)` 
  - `extract_seq_for_side(member, genome, side, anchor, flank)`
  - `process_cluster_side(...)` returning per-member dicts
  - `validate_motif(seq, side, motif)` — central place for TG/CA / none / future
- Unit tests for: motif validation, boundary geometry on +/- strand,
  edge cases (chromosome boundaries, anchor truncation).

### M2 — refine_boundaries.py + refined GFF3 (4–5 days)

- New script `utils/refine_boundaries.py`:
  - Parses GFF3, builds parent → 5'/3' index.
  - MMseqs2 cluster per lineage (reuse from `run_tecap_only.R` /
    `run_parasail_boundary.py`).
  - Per cluster: parasail anchored extension via M1 helpers.
  - Computes per-cluster validation rate.
  - **No MAFFT fallback yet** (stub returns parasail result).
  - Emits Output A (refined GFF3) and per-element TSV.
- New CLI `dante_ltr_refine` (Python wrapper).
- End-to-end test on Alyr: refined GFF3 has expected count of
  refined elements; per-element TSV has correct schema.

### M3 — MAFFT fallback (3–4 days)

- Extract MAFFT change-point logic from `build_ltr_library.R` into a
  helper script `utils/refine_mafft_fallback.R` callable from Python:
  ```
  Rscript utils/refine_mafft_fallback.R --cluster_fa ... --output_tsv ...
  ```
- `refine_boundaries.py` calls it for clusters where parasail
  validation rate < `--mafft_fallback_threshold` and motif=TG/CA.
- Per-element merge logic: prefer parasail when validated; MAFFT
  fallback when parasail's per-element call wasn't validated AND
  MAFFT's was.
- Update per-element TSV to expose both columns and the merge.

### M4 — build_ltr_library.R refined-aware (3 days)

- Add `--refined_gff3` and `--min_validated_members`.
- When `--refined_gff3` given:
  - Read per-element refinement TSV (or parse the GFF3 attributes).
  - Build per-cluster consensus from validated members only.
  - If validated members < `--min_validated_members`, build from all
    members and flag `low_confidence`.
- Emit `low_confidence` flag in the map TSV.
- Backward compat: if `--refined_gff3` not given, the script behaves
  exactly as today.

### M5 — dante_ltr_solo integration (2 days)

- `dante_ltr_solo` accepts `--refined_gff3`.
- Pass through to `build_ltr_library.R` and to coordinate-mapping logic.
- No change to solo-LTR detection algorithm.
- End-to-end test: refined run produces solo-LTR GFF3 of expected size
  on Alyr.

### M6 — configuration polish + docs (2 days)

- Plumb `--boundary_motif` end-to-end.
- For v1, only `TG/CA` and `none` accepted. Validate at CLI parse
  time.
- Future: per-lineage `boundary_motif` column in
  `lineage_domain_order.csv` is a documented stub but not consumed in
  v1 (a TODO in the code).
- Update `README.md` with new command + brief example.

### M7 — tests + benchmarks (2 days)

- Add `tests/refine.sh` exercising the full pipeline on
  `test_data/sample_DANTE_part.fasta` if a refined-LTR-bearing test
  fixture exists, or on Alyr cached data.
- Performance regression test: pipeline on Alyr < 10 min wall.

**Total estimate: ~3 weeks** if no surprises.

## 9. Risks and mitigations

| risk | mitigation |
|---|---|
| Parasail prototype has untested edge cases (chromosome edges, very-short LTRs, chimeras) | Add unit tests in M1 covering: chrom-end clipping, anchor < 10 bp, mid-cluster reverse-complement mismatches |
| MAFFT fallback introduces double-spend on big Pisum-like genomes | Lineage-aware threshold; cap `mafft_fallback_threshold` so fallback fires on at most ~10 % of clusters |
| Refined-GFF3 attribute bloat breaks tools that only read `start`/`end` | Keep attributes additive; don't change column 1–8 |
| Solo-LTR pipeline has hardcoded assumptions about original GFF3 layout | M5 isolates the change behind an `--refined_gff3` switch and tests both paths |
| Lineages without TG/CA give artifactually low validation rates | Document `--boundary_motif none` and how to use it; recommend per-lineage config (M6 stub) for v1.x |
| Per-element refinement adds memory pressure on Pisum-scale | Process clusters in streaming batches; never hold all members in memory at once |

## 10. Open questions

Not blocking M1–M5; revisit before merging.

1. **Snap-window size**: 5 bp default is good on Alyr. Should it scale
   with anchor length or stay constant? — defer; keep 5 as default,
   reassess after running on more genomes.
2. **Per-lineage motif config rollout**: when does this become v1.1?
   Probably after we have observations from multiple genomes where
   global TG/CA causes issues. Add the column to
   `lineage_domain_order.csv` as a stub now (ignored at runtime),
   wire it through in v1.1.

## 11. Validation strategy

After M5 lands:

| stage | dataset | metric |
|---|---|---|
| smoke | `test_data/sample_DANTE_part` | refined GFF3 valid; pipeline runs end-to-end |
| Alyr (g1) | `test_data/g1_dante_ltr.gff3` | per-member TG/CA-OK ≥ 75 % (matches prototype) |
| g2 | `test_data/g2/` | per-member TG/CA-OK ≥ 85 % |
| g3 (Pisum) | `test_data/g3/` (when ceph available) | wall < 4 h; library size in line with previous runs |
| solo-LTR sanity | g1 + g2 | new solo-LTR count within ±20 % of current production |

The "±20 %" range allows for the new pipeline to find both more (better
boundaries → more search hits) and fewer (cleaner library → fewer
spurious hits) solo LTRs; outside that range we'd want to investigate.

## 12. What this plan deliberately defers

- Replacing solo-LTR detection algorithm itself (still BLAST-based).
- Adding non-TG/CA motif support beyond the `none` mode.
- Per-lineage config consumption (kept as stub).
- Telling users to re-run DANTE_LTR with refined boundaries fed back —
  that's a longer-horizon goal mentioned in earlier notes; remains
  separate from this plan.
- HTML report for refinement results — could borrow from
  `utils/boundary_report.R` but is not v1 scope.
