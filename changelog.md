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


