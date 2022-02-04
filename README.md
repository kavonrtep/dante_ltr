# dante_ltr

Tool for identification of complete LTR retrotransposons based on analysis of protein
domains identified by DANTE tool.

## Installation:

```shell
conda create -n dante_ltr -c bioconda -c conda-forge -c petrnovak dante_ltr
```
## Usage

```shell
Usage: ./extract_putative_ltr.R COMMAND [OPTIONS]


Options:
        -g GFF3, --gff3=GFF3
                gff3 with dante results

        -s REFERENCE_SEQUENCE, --reference_sequence=REFERENCE_SEQUENCE
                reference sequence as fasta

        -o OUTPUT, --output=OUTPUT
                output file path and prefix

        -c NUMBER, --cpu=NUMBER
                Number of cpu to use [default 5]

        -h, --help
                Show this help message and exit
```

## Example
```shell
mkdir -p tmp
./extract_putative_ltr.R -g test_data/sample_DANTE.gff3 -s test_data/sample_genome.fasta -o tmp/ltr_annotation
```

## Output files : TODO 