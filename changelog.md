## 0.5.4.0 (2026-08-19)

Chunk-pool memory budget (issue #13).

* New `dante_ltr --max_memory <GB>` flag (alias `--max-memory`) giving
  the memory allocation this run may use, for sizing the per-chunk
  detection pool.  Matches TideCluster 1.20.0's flag name, units and
  resolution chain.
* Without the flag the budget is resolved from, first hit wins:
  `AGENT_MEMORY`; the scheduler variables `PBS_RESC_MEM`,
  `SLURM_MEM_PER_NODE`, `LSB_MAX_MEM_RUSAGE`, `SLURM_MEM_PER_CPU` x
  `SLURM_CPUS_ON_NODE`; the cgroup v2/v1 memory limit, walking
  `/proc/self/cgroup` leaf to root so a limit on an ancestor job scope
  is found; and finally `/proc/meminfo` `MemAvailable` as before.
  `MemAvailable` is not namespaced, so under a batch scheduler or in a
  container it reports the whole node and the pool could be sized many
  times too large and OOM-killed hours into a run.
* The pool-sizing message now names the budget's source, and a warning
  is printed when the budget is host-wide while a scheduler job is
  detected (`PBS_JOBID`, `SLURM_JOB_ID`, `SLURM_JOBID`, `LSB_JOBID`).
* Behaviour on a plain host with no scheduler or container is
  unchanged.  Where a new source applies the only effect is a smaller
  pool; per-chunk results are keyed by index, so the output stays
  byte-identical.
* New test level `./tests.sh mem` (`tests/test_mem_budget.py`), also
  wired into the CI tests workflow.

## 0.5.3.0 (2026-08-10)

Parallel per-chunk LTR detection.

* The large-genome chunk loop in `dante_ltr` ran the independent
  `detect_putative_ltr.R` invocations serially, so the stage used
  ~1 core regardless of `-c`.  Chunks now run in a
  `multiprocessing.Pool` of single-threaded (`-c 1`) children.
* The pool is sized by cores **and** memory: the largest chunk runs
  first as a probe, its peak RSS is read via
  `getrusage(RUSAGE_CHILDREN)`, and the pool width is capped at
  `0.8 * MemAvailable / peak` so a large `-c` cannot OOM the host.
  (The budget term is corrected in 0.5.4.0.)
* Chunks are dispatched largest-first to drop the straggler tail,
  workers ignore `SIGINT` so Ctrl-C is handled once in the parent, and
  the fd budget reserves descriptors for the concurrent children.
* Output is unchanged — every downstream step walks the chunks in
  index order, so completion order is never observed.

## 0.5.2.0 (2026-07-21)

Deterministic LTR-RT library clustering.

* `mmseqs easy-cluster` is order-sensitive: the same sequence set in a
  different order elects different representatives and a different
  cluster count.  `TE_all.fasta` arrived in `DANTE_LTR.gff3` row order,
  which is not canonical on the multi-chunk path, so the library — and
  the downstream RepeatMasker annotation — varied run-to-run and across
  machines (35 vs 36 representatives on a 223-sequence fixture).
* `utils/mmseq_clustering.R` now sorts `TE_all` into a canonical order
  by sequence content (then name, to break ties) before the partition
  step that feeds mmseqs, making the library a deterministic function
  of the input *set* regardless of record order, and robust to upstream
  coordinate jitter.

## 0.5.1.1 (2026-07-18)

Large-genome performance in the Python wrapper.

* Removed the quadratic scans from post-R GFF3 assembly, which
  dominated runtime on large genomes (hours on a 90 Gbp assembly).
  The matching table is indexed once by header instead of being
  rescanned per feature line; `split_fasta_to_chunks()` uses the same
  index in place of its O(contigs²) scan.
* Dropped per-line `Gff3Feature` construction from the coordinate
  recalculation and filtering passes (4 construction sites → 1).
  `normalize_gff3_attributes()` reproduces the attribute normalisation
  the old dict rebuild performed implicitly, so output is unchanged.
* `split_fasta_to_chunks()` streams the reference a record at a time
  instead of loading it into a dict of Python strings.  Peak RSS is now
  flat at ~14 MB for 7.7 / 86 / 244 Mbp inputs (was 35 / 248 / 676 MB),
  i.e. interpreter baseline, independent of genome size.
* Tests cover the `dante_ltr_solo` chunked path.

## 0.5.1.0 (2026-07-14)

Refinement v2.1, the solo pipeline on refined GFF3, and the
"Too many open files" fix.

* **Boundary refinement v2.1**
  * New tier-3 per-side MSA rescue (`--msa_rescue`, on by default):
    MAFFT runs on every qualifying cluster, and per-side rescue applies
    to any side still unresolved after the earlier tiers.  The TSD-loss
    revert rule gates these proposals like the others.
  * Confidence split into two axes.  `Refinement_Status`
    (`not_evaluated` / `unresolved` / `confirmed` / `refined`) answers
    "what happened to this side?"; `Refinement_Confidence`
    (`dual` / `divergent` / `inner_only` / `mafft` / `unrefined`)
    records the source and quality.
  * MAFFT snap window widened 5 → 20 bp: the change-point detector
    fires 5-15 bp inside the conserved region, so the old window failed
    motif validation on ~10 % of MSA calls that were within ±20 bp of a
    real TG/CA.
  * Fixed MAFFT `motif_ok=FALSE` on genuinely matching motifs —
    `as.character(DNAStringSet)` returns a *named* character vector and
    `identical()` compares attributes, so every MAFFT-derived boundary
    was rejected.  Restores `motif_ok=TRUE` on ~94 % of MAFFT calls.
  * `Cluster_ID` / `Cluster_Size` are now emitted for **all** clustered
    members, with `Refinement_Status=not_evaluated` on those below
    `--min_cluster_size`, making refinement the single source of
    clustering truth.
* **Solo LTR pipeline**
  * A refined GFF3 is now the universal input; an unrefined DANTE_LTR
    GFF3 triggers `dante_ltr_refine` automatically, written to
    `<output_dir>/refined/`.  This removes the `--refined_gff3` and
    `--min_validated_members` flags added in 0.4.1.0.
  * `utils/build_ltr_library.R` drops its own MMseqs2 clustering and
    its cluster-wide TG/CA + TSD validation, grouping by `Cluster_ID`
    and ranking within a cluster by `Refinement_Status` instead.
  * `LibraryConfidence` (`validated` / `mixed` / `unrefined`) is
    propagated onto every emitted `solo_LTR` feature.
  * `SL_noTSD` representatives showing a positive `UTR5_junction`,
    `PPT_junction` or `PBS_check` signal are reclassified as likely
    fragments of complete elements and partitioned into
    `solo_ltr_te_fragments.gff3` with a parallel statistics CSV.  On the
    *A. lyrata* test run this moved 2100 features out of a
    3587-feature `solo_ltr.gff3`.
* **Fixed "Too many open files" on large genomes** (#12).  The chunk
  split opened one handle per chunk at once — ~1800 for a 90 Gbp genome
  at `-S 50000000`, against a default soft limit of 1024.
  `_limit_chunks_to_fd_budget()` raises the soft limit toward the hard
  limit and only then reduces the chunk count, so the run degrades
  gracefully instead of aborting.  Adds `./tests.sh fd`.
* **Docs**: README restructured — feature sections compressed, CLI
  reference moved to the bottom.
* **CI**: release workflow installs runtime deps from `conda-deps.txt`
  and builds with `conda build` from the base env (drops boa).

## 0.5.0.0 (2026-05-06)

Refinement v2 — inner-primary policy.  (Not tagged; first published in
0.5.1.0.)

* Two parasail passes per query × side: an outer-edge pool (same role)
  and an inner-edge pool (opposite role used at the inner edge).  The
  inner pool's flank asymmetry sharpens the boundary, and same-role
  drift is no longer the primary signal.
* Inner-primary policy: take the inner coordinate when its TG/CA motif
  validates, otherwise keep DANTE_LTR's original.
* Per-element TSD recheck at the proposed coordinates, with a revert
  rule — if applying a proposal destroys the originally detected TSD,
  every changed side on that element reverts to the original.
* New `utils/tsd_redetect.py`, a pure-Python port of the
  `evaluate_ltr()` TSD logic from `utils/ltr_utils.R`, validated
  bit-for-bit against the R implementation.
* `utils/parasail_boundary.py` gains `extract_pool_for_side()`,
  `aggregate_extension()` and `project_corrected_g()` so a caller can
  pass an arbitrary anchor pool.
* Design: `docs/refine_v2_analysis.md`.

## 0.4.1.0 (2026-05-06)

Per-element LTR boundary refinement.  (Not tagged; first published in
0.5.1.0.)

* New `dante_ltr_refine` command — parasail anchored extension across
  MMseqs2-clustered members, with optional MAFFT change-point fallback
  for low-confidence clusters.  Outputs a refined GFF3 with the same
  schema plus `Original_Start` / `Original_End` /
  `Refinement_Method` / `Refinement_Confidence` /
  `Cluster_ID` / `Cluster_Size` / `TG_OK` / `CA_OK` attributes, plus
  a per-element TSV, cluster manifest TSV and run JSON.
* `dante_ltr_solo` accepts `--refined_gff3 PATH` and
  `--min_validated_members N` (default 4).  When supplied, the LTR
  reference library is built from validated members only;
  low-confidence clusters fall back to all members and are flagged in
  the library map TSV.
* `utils/build_ltr_library.R` accepts `--refined_gff3` (skips the
  per-cluster change-point step in this mode and builds the
  consensus by majority vote over a MAFFT alignment of the validated
  bodies) and `--min_validated_members`.
* New library module `utils/parasail_boundary.py` (importable) with
  the parasail anchored-extension primitives, motif validation and
  per-side / per-strand boundary geometry.
* See `docs/parasail_boundary_refinement_plan.md` for design and
  validation strategy.

Also in this release:

* **Fallback classification fix**: `demote_classification()` returned
  early on level-4 input such as
  `Class_I|LTR|Ty3/gypsy|chromovirus`, so those features were left
  unchanged under `coarse2` and never reached the loose `Ty3/gypsy`
  row.  Rewritten with branch-aware semantics (108 vs 83 TEs on the
  `fallback_medium` fixture).  `--fallback_mode` is now documented in
  the README, flagged experimental.
* **`utils/build_ltr_library.R` boundary work**: tighter change-point
  scan (`conservation_window` 5 → 7, new `high_sustain` requiring two
  consecutive columns above threshold); `min_cluster_size` default
  3 → 6; a wide-flank retry for clusters whose entire flank is
  conserved, capped proportionally at
  `min(--wide_flank, max(500, median_annot_body_len))` with
  `--wide_flank` default 1000; report-only TSD and TG/CA validation of
  each correction (four new QC TSV columns, measured at both the
  annotated and corrected column pair); and a self-contained HTML
  boundary report (`utils/boundary_report.R`, eight panels, PNGs
  inlined as base64, no new R dependencies).
* **Packaging**: `parasail-python` and `pyfaidx` added as run deps and
  `dante_ltr_refine` exposed by the recipe; `requirements.txt` renamed
  to `conda-deps.txt` so GitHub's Dependabot pip detector stops trying
  to parse conda pin syntax.

## 0.4.0.13 (2026-04-20)

Fallback classification mode for distantly-related species.

* New `dante_ltr --fallback_mode {none,coarse3,coarse2}` flag. When
  enabled, LTR `Final_Classification` values are demoted to a coarser
  taxonomic depth before matching, and a correspondingly relaxed
  built-in constraints table is used.
  * `coarse3` — classify into `Ty1/copia`, `Ty3/gypsy/chromovirus`,
    `Ty3/gypsy/non-chromovirus`.
  * `coarse2` — classify into `Ty1/copia`, `Ty3/gypsy` only.
* Every demoted feature carries an `Original_Classification`
  attribute so nothing is lost downstream.
* Always-on **recommendation heuristic**: the Python wrapper logs
  one line summarising LTR-domain resolution (filtered by
  `Relat_Length >= min_relative_length`) and, when
  `--fallback_mode=none`, suggests a fallback mode when the input
  has a lot of under-resolved calls.  No behaviour change beyond
  the log line for well-classified inputs (`sample_DANTE*`: 1 %
  superfamily-only, no recommendation; Drapa subset: 68 %, suggests
  `coarse2`).
* **Data-driven aRH inclusion** for `coarse3`. After demotion, the
  wrapper spatially clusters non-chromovirus LTR domains; when aRH
  appears in ≥ 20 % of clusters, the non-chromovirus row's
  `Domains order` is rewritten to `GAG PROT RT RH aRH INT` in a
  run-time-only copy of the constraints table (written to
  `<out_dir>/_fallback_inputs/constraints_coarse3_arh.csv` for
  inspection).  The shipped constraints tables on disk are never
  modified.
* Two new shipped constraints tables:
  `databases/lineage_domain_order_coarse3.csv` (3 rows) and
  `databases/lineage_domain_order_coarse2.csv` (2 rows).
* New test data `tests/data/fallback_medium/` — 4 MB single-contig
  slice of a Drapa assembly (distantly related, no aRH). New
  `tests/fallback.sh` exercises both the aRH-absent (Drapa) and
  aRH-present (`sample_DANTE_part`) arms end-to-end; wired into
  `tests.sh`, `./tests.sh all`, and the GitHub Actions workflow.

Non-fallback runs (`--fallback_mode=none`, default) on the existing
test data produce the same TE counts as 0.4.0.12 — only a new
informational log line differs.

Design: `docs/fallback_classification_design.md` (draft with user
annotations), `docs/fallback_classification_spec.md` (locked v1
spec), `docs/fallback_classification_implementation_plan.md`.

## 0.4.0.12 (2026-04-18)

Release-workflow fix: pin `setuptools<81` in the release env.

0.4.0.11 built the correct package (`dante_ltr-0.4.0.11-0.tar.bz2`)
but the upload step crashed with
`ModuleNotFoundError: No module named 'pkg_resources'`.
`anaconda-client`'s `binstar_client` still imports `pkg_resources`,
which `setuptools 81.0.0` removed from the default install.

Verified locally: with `setuptools 82.0.1` the `anaconda` CLI fails
to load; with `setuptools 80.10.2` (the last pre-81 release) the CLI
loads cleanly and `anaconda upload --help` prints usage.

Upstream already flags the workaround:
`"pkg_resources is deprecated ... pin to Setuptools<81"`.
The pin stays until `anaconda-client` migrates to
`packaging.version`.

## 0.4.0.11 (2026-04-18)

Another release-workflow fix.  conda-build's `load_file_data`
accepts only JSON/YAML/TOML and raised
`ValueError: Unknown file format: py` during the second rendering
pass, after the first pass had already fallen back to `0.0.0.dev`.

`meta.yaml` now uses `load_file_regex` to extract `__version__` from
the copied `version.py` via a regex pattern, which is the
general-purpose text loader in conda-build's jinja context.
Verified locally that the match yields `0.4.0.11`.

## 0.4.0.10 (2026-04-18)

Further release-workflow fixes. No user-visible functional change.

* `conda/dante_ltr/meta.yaml` is back to reading the version from a
  `version.py` file in the recipe directory via `load_file_data`. The
  env-var-based approach in 0.4.0.9 did not work: conda-build's jinja
  context didn't propagate `PKG_VERSION` from the workflow step's
  `env:` block during rendering, so the build was still producing
  `dante_ltr-None-0.tar.bz2`.
* `.github/workflows/conda-release.yml` now `cp`s the root
  `version.py` into `conda/dante_ltr/` before `conda mambabuild`.
  The copy is `.gitignore`d (canonical version.py stays at the root).
* Release workflow installs `setuptools` alongside `anaconda-client`
  to restore `pkg_resources`, which the client still imports. On
  Python 3.12+ environments without setuptools this manifests as
  `ModuleNotFoundError: No module named 'pkg_resources'` during
  upload.

## 0.4.0.9 (2026-04-18)

Release-workflow fixes for 0.4.0.8. Same user-visible functionality.

* `conda/dante_ltr/meta.yaml` now reads the version from the
  `PKG_VERSION` environment variable instead of
  `load_file_data('../../version.py')`. conda-build copies the recipe
  to a temp workdir before rendering, so relative paths *above* the
  recipe don't resolve and the build was shipping a package named
  `dante_ltr-None-0`.
* `.github/workflows/conda-release.yml` passes `$GITHUB_REF_NAME` to
  the build step as `PKG_VERSION`, matching the tag-as-version sanity
  check earlier in the workflow.
* `tests/long.sh` defaults to the committed `sample_genome.fasta`; the
  g1 fallback warning is gone. Override via `LONG_FASTA` / `LONG_DANTE`
  env vars for local runs on bigger datasets.

## 0.4.0.8 (2026-04-17)

Infrastructure release — first version published through the new
GitHub-Actions conda release workflow.

* **CI** (`.github/workflows/`):
  * `tests.yml` — smoke + short tests on every push and PR.
  * `conda-release.yml` — on tag push, runs the long test as a gate,
    builds the conda package via `conda mambabuild`, and uploads to the
    `petrnovak` anaconda channel. Requires the `ANACONDA_API_TOKEN`
    repo secret.
* **Tests** (`tests/`):
  * `smoke.sh` — < 30 s, exercises every CLI and a tiny pipeline on a
    40 kb window (`tests/data/smoke/`).
  * `short.sh` — ~ 1 min, full pipeline on the existing
    `test_data/sample_genome_part.fasta`.
  * `long.sh` — ~ 10-30 min, full pipeline on `test_data/g1*` (falls
    back to `sample_genome.fasta` when `g1` is not present).
  * Root `tests.sh` becomes a dispatcher; `./tests.sh <N>` still runs
    the long test with `N` CPUs for backwards compatibility.
* **In-repo conda recipe** (`conda/dante_ltr/`): `meta.yaml` reads the
  version from `version.py` via jinja, sources from the local tree
  (`path: ../..`, no sha256 dance); `build.sh` copies the tree and
  symlinks every CLI into `$PREFIX/bin`. Adds the `mafft` runtime
  dependency (missing from the previous external recipe) and the
  `dante_ltr_solo` / `dante_ltr_gff3_to_canonical` / `dante_reclassify`
  executables.
* **Docs** moved to `docs/`, keeping the repo root clean (README and
  changelog stay at the root).
* **Bug fix**: `select_solo_representatives.R` now emits an empty
  statistics CSV when the input GFF3 has zero `solo_LTR` features,
  instead of early-exiting without the file.

## 0.4.0.7 (2026-04-17)

Joint 5'+3' LTR MSA for true boundary-aware consensus.

* **`utils/build_ltr_library.R`**:
  * MMseqs2 clustering still runs on 5'LTR bodies, but each cluster's
    MAFFT input now also includes the sibling 3'LTRs (resolved by the
    shared `Parent` TE id).
  * Per-column conservation is now computed separately for the 5'LTR
    subset and the 3'LTR subset. The 5' boundary is a change-point on
    the 5'LTR-subset profile (random genomic 5' flank); the 3' boundary
    is a change-point on the 3'LTR-subset profile (random genomic 3'
    flank). Previously the 3' boundary fell back to the median
    annotated column.
  * Consensus is the column-wise majority across the full joint matrix
    between the two corrected boundaries — both subsets contribute,
    reinforcing the body consensus.
  * New helper `detect_3p_change_point()` mirrors `detect_5p_change_point()`.
  * `mafft_boundary_consensus()` gains a `roles` parameter; backwards-
    compatible when `NULL` (falls back to median-annotated 3' end).
* **Boundary QC TSV** gains `n_5ltr`, `n_3ltr`, `corrected_3_col`,
  `shift_3`.

Impact on reference test (`test_data/g1_dante_ltr.gff3`):
* 888 consensuses in both before and after; median length unchanged
  (389 → 390 bp).
* `shift_5` distribution unchanged (same data, same algorithm).
* `shift_3` now non-zero in 122/180 multi-member clusters (median +4 bp,
  max +49 bp) — all extensions rightward into the annotated 3' flank,
  consistent with DANTE_LTR annotations being slightly conservative on
  the 3' side.
* MAFFT wall time ~2× because the input per cluster now carries both
  LTRs (library build 186 s → 299 s end-to-end).

## 0.4.0.6 (2026-04-17)

`dante_ltr_solo` pipeline overhaul — library build, TSD detection,
overlap handling, and performance.

* **Library build (`utils/build_ltr_library.R`)**:
  * 5'LTR-only flank-aware consensus. LTRs are extracted with `±flank` bp
    (default 50) and aligned with MAFFT; the 5' boundary is refined by a
    change-point scan on the per-column conservation profile. The 3'
    boundary uses the median annotated position across the cluster.
  * 3'LTRs are still used for the PPT tag database but no longer enter the
    library consensus.
  * New per-cluster QC TSV `*_LTR_library_boundary_qc.tsv` with annotated
    vs corrected boundary columns, shift magnitudes, and body lengths.
* **TSD detection (`check_tsd` in `utils/solo_ltr_utils.R`)**:
  * Scans `(mode-1):(mode+1)` around the per-lineage modal TSD length
    (clipped to 3-8); falls back to 4-6 bp when no lineage entry.
  * 1-mismatch tolerance relaxed to length ≥ 4 (was ≥ 5).
  * Two-pass search: exact match (longest first) always wins over 1-mm.
* **LibraryID on every `solo_LTR`**: every hit carries the library
  consensus id (`LibraryID=LTR_XXXXXX`) that produced it.
* **Representative-per-locus output**: overlapping hits are now collapsed
  with a ≥ 50 % reciprocal-overlap graph; two GFF3 files are produced:
  * `solo_ltr_raw.gff3` — all hits (deduped by coordinate)
  * `solo_ltr.gff3` — one representative per locus (SL preferred over
    SL_noTSD, then longest, tie-break by identity/LibraryID)
  * Representatives carry `ClusterSize`, `SupportingHits`, plus optional
    `boundary_uncertain=true` (SL nested inside a longer SL_noTSD) and
    `class_conflict=true` flags. New script
    `utils/select_solo_representatives.R`.
* **Performance**:
  * Vectorised `make_solo_ltr_gff3()` (~260× faster on large chunks).
  * Batched junction/PBS BLAST — one blastn invocation per database
    instead of three per hit (~38× faster on that stage).
  * Vectorised MSA helpers (conservation, sliding mean, consensus,
    position mapping) and union-find for cluster grouping.
  * Per-stage wall-clock profiling printed by every R script.
  * End-to-end `detect_solo_ltr` ~4.6× faster (835 s → 180 s on the
    reference test).

## 0.4.0.1 (2024-09-20)

*  fix in installation instruction 
*  updated to REXdb Viridiplantae 4.0 

## 0.3.5.1 (2024-04-18)

*  Improved HTML report 


## 0.3.5.0 (2023-11-28)
 
* Option to ignore ambiguous domains added 
* Report is now in html format 

## 0.3.4.0 (2023-11-22)

*  bug fix in TSD reporting. For TE on minus strand, TSD sequence was reported as 
   sequence on plus strand - this is fixed 
*  possibility to adjust te_detection constrains added, TE constrains updated
   documentation updated 
*  reporting flanking sequences around TE documentation updated 


## 0.2.3.8 (2023-10-23)

*  bugfix in urldecoding, summary script updated 

## 0.2.3.7 (2023-10-13)

*  dante_ltr_summary export csv table with basic stats 
*  bugfix - overlapping domain from DANTE handling added 

## 0.2.3.6 (2023-09-25)

*  docs updated 
*  bugfix in dante_ltr_summary (error when no element in the lineage was complete) 
*  
## 0.2.3.5 (2023-08-31)

*  bug fix in path to ltr_utils.R 


## 0.2.3.4 (2023-08-16)

*  bugfix - no domain contigs handling 


## 0.2.3.3 (2023-08-10)

*  fix in gff3 output - source attribut changed from ../dante_ltr to dante_ltr 

## 0.2.3.2 (2023-07-25)

*  bugfix in TSD identification 

## 0.2.3.1 (2023-07-25)

*  help updated 

## 0.2.3.0 (2023-07-25)

*  bug fix  #4 
*  new script added - extract elements as fasta from gff 
*  version writen to gff output #5 
*  dante_ltr_to_library added 

## 0.2.2.4 (2023-05-05)

*  docs updated 
*  bug fix in danter_ltr summary 

## 0.2.2.3 (2023-05-03)

*  version printing added  
*  readme updated 
*  summary script added

## 0.2.2.2 (2023-03-24)

*  bugfix - error when no TE detected 
*  bugfix #3 


## 0.2.2.1 (2023-01-31)

*  readme updated 
*  readme updated, some print statement removed 
*  save complete environment if debug is on 


## 0.2.2.0 (2023-01-23)

*  chdcr bugfix 
*  gitpod conf added 
*  bugfix trna database path + chdcr 
*  gitpod updated 
*  bug fix in path for trna_db_hemi 


## 0.2.1.0 (2023-01-03)

*  readme updated 
*  Ndomains attribute added to transposable_element 
*  Ndomains attribute added to transposable_element - bugfix 
*  allow for one mismatch in TSD 
*  bug fix in calculation of domain distances, debug option added 
*  bug fix in lineage ltr constrains 
*  improved TG/CA detection mismatches are not countexd into TG/CA distance from alignment end 
*  adjusted threshold for dante filtering, ltr distance constrains adjusted 
*  improved removing of overlaping partial te, corrected reporting of LTR coordinated (some LTR did not have TG/CA) 
*  database of tRNA updated - 1/2tRNA added 
*  evalue for PBS added, some bug fixes 
*  evalue changed to PBS_evalue 


## 0.1.9.0 (2022-11-30)

*  change log added 
*  tandem repeat masking - using TideHunter output gff is not sorted 
*  validation of input added 


## 0.1.8.0 (2022-08-10)

*  dante_ltr now runs on chunks if input is large, long contigs split to chunks too more test files added 


## 0.1.7.0 (2022-07-21)

*  version in galaxy files updated 
*  bugfix - correction in get_ranges functions 
*  bugfix -handling sequences with no DANTE domains 
*  dante_ltr now runs on chunks if input is large 


## 0.1.6.1 (2022-06-28)

*  improved detection of incomplete elements, better sensitivity statistical model for individual lineages added 
*  documentation update 
*  documentation update + figure 
*  filtering of DANTE is part of preprocessing of DANTE GFF3 this allow to use full DANTE gff -> better sensitivity 


## 0.1.6.0 (2022-06-24)

*  improved detection of incomplete elements, better sensitivity statistical model for individual lineages added 


## 0.1.5.4 (2022-05-31)

*  new option added - to tolerate missing domains 


## 0.1.5.3 (2022-05-19)

*  bugfix - adjusted orge and retand max ltr distance 
*  tool description updated added common ID to fasta libraries output 


## 0.1.5.2 (2022-05-16)

*  bugfix - adjusted athila max ltr distance 


## 0.1.5.1 (2022-05-09)

*  bugfix - handle duplicated ID on input 


## 0.1.5 (2022-05-03)

*  version in tool xml updated, typo in xml corrected 
*  full TE output added, summary pdf output added 


## 0.1.4 (2022-04-12)

*  galaxy tool definition added 
*  documentation updated 
*  bug fix - typo in requirements 
*  search of ltr boundaries adjusted - more sensitive 
*  version in tool xml updated 


## 0.1.3.1 (2022-02-04)

*  singularity definition corrected 
*  singularity definition corrected 
*  singularity definition corrected 
*  bug fix - ltr_utils added to git 


## 0.1.3 (2022-02-04)

*  dokumentation updated - installation instruction using conda 
*  R function moved to separate file 
*  bug fixes, script for annotation cleaning added 


## 0.1.2 (2022-01-05)

*  bugfix - path to alternative configuration corrected 


