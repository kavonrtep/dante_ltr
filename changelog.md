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


