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


