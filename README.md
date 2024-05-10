# DANTE_LTR


[![Anaconda-Server Badge](https://anaconda.org/petrnovak/dante_ltr/badges/version.svg)](https://anaconda.org/petrnovak/dante_ltr)  [![DOI](https://zenodo.org/badge/439021837.svg)](https://zenodo.org/badge/latestdoi/439021837)

Tool for identifying complete LTR retrotransposons based on analysis of protein domains identified with the [DANTE tool](https://github.com/kavonrtep/dante). Both DANTE and DANTE_LTR are available on [Galaxy server](ttps://repeatexplorer-elixir.cerit-sc.cz/).

## Preprint
DANTE and DANTE_LTR: computational pipelines implementing lineage-centered annotation of LTR-retrotransposons in plant genomes, Petr Novak, Nina Hostakova, Pavel Neumann, Jiri Macas
bioRxiv 2024.04.17.589915; doi: https://doi.org/10.1101/2024.04.17.589915

## Principle of DANTE _LTR
Complete retrotransposons are identified as clusters of protein domains recognized by the DANTE tool. The domains in the clusters must be assigned to a single retrotransposon lineage by DANTE. In addition, the orientation and order of the protein domains, as well as the distances between them, must conform to the characteristics of elements from REXdb database [Neumann et al. (2019)](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0144-1). 
In the next step, the 5' and 3' regions of the putative retrotransposon  are examined for the presence of 5' and 3' long terminal repeats. If 5'- and 3'-long terminal repeats are detected, detection of target site duplication (TSD) and primer binding site (PSB) is performed. The detected LTR retrotranspsons are classified into 5 categories:
- Elements with protein domains, 5'LTR, 3'LTR, TSD and PBS - rank **DLTP**.
- Elements with protein domains, 5'LTR, 3'LTR, and PBS (TSD was not found) Rank **DLP**
- Elements with protein domains, 5' LTR, 3'LTR, TSD (PBS was not found) - rank **DTL**
- Elements with protein domains, 5'LTR and 3'LTR (PBS and TDS were not found) - rank **DL**
- Elements as clusters of protein domains with the same classification, no LTRs - rank **D**.

![dante_ltr_workflow.png](dante_ltr_workflow.png)

## Availability
DANTE_LTR and DANTE are available on [Galaxy server](https://repeatexplorer-elixir.cerit-sc.cz/) or can be installed using conda package manager.

## Installation:
[![Anaconda-Server Badge](
https://anaconda.org/petrnovak/dante_ltr/badges/version.svg)](https://anaconda.org/petrnovak/dante_ltr)

```shell
conda create -n dante_ltr -c bioconda -c conda-forge -c petrnovak dante_ltr
```

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/kavonrtep/dante_ltr)

## Input data
One input is a reference sequence in fasta format. The second input is an annotation of the reference genome using the [DANTE tool](https://github.com/kavonrtep/dante). For better results, use the unfiltered full output of the DANTE pipeline.


## Usage

### Detection of complete LTR retrotransposons

```
usage: dante_ltr [-h] -g GFF3 -s REFERENCE_SEQUENCE -o OUTPUT [-c CPU]
                 [-M MAX_MISSING_DOMAINS] [-L MIN_RELATIVE_LENGTH] [-S MAX_CHUNK_SIZE]
                 [-v] [--te_constrains TE_CONSTRAINS] [--no_ambiguous_domains]

        Tool for identifying complete LTR retrotransposons based on 
        analysis of protein domains identified with the DANTE tool
        

options:
  -h, --help            show this help message and exit
  -g GFF3, --gff3 GFF3  gff3 file with full output from Domain Based Annotation of Transposable Elements (DANTE)
  -s REFERENCE_SEQUENCE, --reference_sequence REFERENCE_SEQUENCE
                        reference sequence as fasta file
  -o OUTPUT, --output OUTPUT
                        output file path and prefix
  -c CPU, --cpu CPU     number of CPUs
  -M MAX_MISSING_DOMAINS, --max_missing_domains MAX_MISSING_DOMAINS
  -L MIN_RELATIVE_LENGTH, --min_relative_length MIN_RELATIVE_LENGTH
                        Minimum relative length of protein domain to be considered for retrostransposon detection
  -S MAX_CHUNK_SIZE, --max_chunk_size MAX_CHUNK_SIZE
                        
                                If size of reference sequence is greater than this value, reference is '
                                'analyzed in chunks of this size. default is 100000000 '
                                'Setting this value too small  will slow down the analysis
                                
  -v, --version         show program's version number and exit
  --te_constrains TE_CONSTRAINS
                        csv table specifying TE constraints for LTR search, template for this table 
                        can be found in https://github.com/kavonrtep/dante_ltr/blob/main/databases/lineage_domain_order.csv
  --no_ambiguous_domains
                        Remove ambiguous domains from analysis
```

#### Example:

```bash
mkdir -p tmp
./dante_ltr -g test_data/sample_DANTE.gff3 -s test_data/sample_genome.fasta -o tmp/ltr_annotation
```

####  Files in the output of `extract_putative_ltr.R`:

- `prefix.gff3` - annotation of all identified elements
- `prefix_D.fasta` - partial elements with protein **d**omains
- `prefix_DL.fasta` - elements with protein **d**omains and **L**TR
- `prefix_DLTP.fasta` - elements with **d**omains, **L**TR, **T**SD and **P**BS
- `prefix_DLP.fasta` - elements with **d**omains, **L**TR and **P**BS
- `prefix_DLT.fasta` - elements with **d**omains, **L**TR, **T**SD 
- `prefix_statistics.csv` - number of elements in individual categories  
- `prefix_summary.html` - graphical summary of the results
- 


### Making library of LTR RTs for  RepeatMasker

If you want to annotate LTR RT elements with custom library using similarity based approach, you can use 
`dante_ltr_to_library` script wich will create non-redundant library which is 
formatted for RepeatMasker:

``` 
usage: dante_ltr_to_library [-h] -g GFF3 -s REFERENCE_SEQUENCE -o OUTPUT_DIR [-m MIN_COVERAGE] [-c CPU]

Creation of repeat library from dante_ltr output. Extract sequences based on gff3 inpute and reference fasta file. Run mmseqs2 clustering to cluster similar sequences to reduce library size. Exclude
clusters which have conflicting annotations and coverage below specified threshold.

options:
  -h, --help            show this help message and exit
  -g GFF3, --gff3 GFF3  gff3 file
  -s REFERENCE_SEQUENCE, --reference_sequence REFERENCE_SEQUENCE
                        fasta file
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory
  -m MIN_COVERAGE, --min_coverage MIN_COVERAGE
                        Minimum coverage of cluster to be included in repeat library (default: 3)
  -c CPU, --cpu CPU     Number of cpus to use

```


### Example of complete workflow of annotation of LTR RTs using DANTE and DANTE_LTR

#### Installation of both DANTE and DANTE_LTR using conda into single environment:
```shell
conda create dante_ltr -c bioconda -c conda-forge -c petrnovak dante_ltr dante
conda activate dante_ltr
```
#### Download example data:

```shell
wget https://raw.githubusercontent.com/kavonrtep/dante_ltr/main/test_data/sample_genome.fasta
```

##### Run DANTE on sample genome using 10 cpus:
```shell
dante -q sample_genome.fasta -o DANTE_output.gff3 -c 10
```
Output will contain annotation of individual protein domains identified by DANTE stored in GFF3 file named `DANTE_output.gff3`. Check DANTE documentation for more details (https://github.com/kavonrtep/dante)  

#### Identify complete LTR retrotransposons  from DANTE ouput using DANTE_LTR
```shell
dante_ltr -g DANTE_output.gff3 -s sample_genome.fasta -o DANTE_LTR_annotation -M 1
```
Option `-M 1` will allow one missing domain in the complete LTR retrotransposon. 

Output files will include: 
- `DANTE_LTR_annotation.gff3` - Annotation of all identified elements
- `DANTE_LTR_annotation_statistics.csv` - number of elements in individual categories
- `DANTE_LTR_annotation_summary.html` - graphical summary of the results

#### Create library of LTR RTs for similarity based annotation
```shell
dante_ltr_to_library -g DANTE_LTR_annotation.gff3 -s sample_genome.fasta -o LTR_RT_library.fasta
```
This step will create non-redundant library of LTR-RT sequences suitable for similarity based annotation using RepeatMasker.


### GFF3 DANTE_LTR output specification
Types of features in GFF3:
- **target_site_duplication**: This feature represents the direct repeats of host DNA 
produced at the insertion site of a transposable element.
- **transposable_element**: This is the main feature representing the full extent of a 
  transposable element within the genome.
- **long_terminal_repeat** (LTR): These are the repetitive sequences found at both ends of 
  retrotransposons (a type of transposable element).
- **protein_domain**: This feature indicates a specific domain within a protein that is 
  part of a transposable element. This corresponds to domains identified by DANTE.
- **primer_binding_site**: This feature represents the site where a primer binds to 
  initiate reverse transcription, typically found in retroviruses and retrotransposons.

Attributes of features in GFF3:
- **Rank**: Rank of the elements (D, DL, DLT, DLP, DLTP) as described above
- **Parent**: Indicates the parent feature of the current feature, here a 
  transposable element ID.
- **Ndomains**: The number of protein domains found within a transposable element.
- **ID**: A unique identifier for the feature.
- **LTR_Identity**: The percentage identity of the LTR sequences.
- **LTR5_length** and **LTR3_length**: The lengths of the 5' and 3' LTRs, respectively.
- **TSD (Target Site Duplication)**: The sequence of the target site duplication.
- **Final_Classification**: A hierarchical classification of the transposable element 
  based on REXdb classification system
- **Name**: The attibute is part of DANTE output and correspod to name of protein 
  domain (RT. RH, PROT, ...)
- **trna_id**: The identifier for the tRNA related to the primer binding site.
- **PBS_evalue**: The E-value associated with the primer binding site.
- **Best_Hit**: Information about the best match of the protein domain to a known database entry.
- **Best_Hit_DB_Pos**: Position of the best hit within the database.
- **DB_Seq**: The database sequence that corresponds to the best hit.
- **Region_Seq**: The sequence of the region in the query that corresponds to the best hit.
- **Query_Seq**: The sequence of the query used to find the best hit.
- **Identity**: The percentage identity of the best hit match.
- **Similarity**: The similarity score of the best hit match.
- **Relat_Length**: The relative length of the match compared to the database sequence.
- **Relat_Interruptions**: Indicates the relative number of interruptions in the 
  domain sequence.Interuption could be either stop codon or frameshift.
- **Hit_to_DB_Length**: The length of the hit compared to the database sequence length.
