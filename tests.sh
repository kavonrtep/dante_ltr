#!/bin/bash
# set to stop in case of error
set -e
# first argument is cpu number
NCPU_TO_USE=$1
eval "$(conda shell.bash hook)"
conda activate dante_ltr
echo "Running tests 1, detection of LTRs"
./extract_putative_ltr.R -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output1 -c $NCPU_TO_USE


cat tmp/test_output1_statistics.csv

echo "Running tests 2, filtering gff"
./clean_ltr.R -g tmp/test_output1.gff3 -s test_data/sample_genome.fasta \
-o tmp/test_output2 -c $NCPU_TO_USE

echo "Running tests 3, detection of LTRs, allow missing domains"
./extract_putative_ltr.R -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output3 -c $NCPU_TO_USE -M 2


cat tmp/test_output3_statistics.csv

echo "Running tests 4, filtering gff"
./clean_ltr.R -g tmp/test_output3.gff3 -s test_data/sample_genome.fasta \
-o tmp/test_output4 -c $NCPU_TO_USE


