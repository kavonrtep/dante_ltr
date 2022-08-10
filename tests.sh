#!/bin/bash
# set to stop in case of error
set -e
# first argument is cpu number
NCPU_TO_USE=$1

# test if NCPU_TO_USE is set:
if [ -z "${NCPU_TO_USE}" ]; then
    echo "NCPU_TO_USE is not set, using 10"
    NCPU_TO_USE=10
fi

eval "$(conda shell.bash hook)"
conda activate dante_ltr
echo "Running tests 1, detection of LTRs"
./detect_putative_ltr.R -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output1 -c $NCPU_TO_USE


cat tmp/test_output1_statistics.csv

echo "Running tests 2, filtering gff"
./clean_ltr.R -g tmp/test_output1.gff3 -s test_data/sample_genome.fasta \
-o tmp/test_output2 -c $NCPU_TO_USE

echo "Running tests 3, detection of LTRs, allow missing domains"
./detect_putative_ltr.R -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output3 -c $NCPU_TO_USE -M 2


cat tmp/test_output3_statistics.csv

echo "Running tests 4, filtering gff"
./clean_ltr.R -g tmp/test_output3.gff3 -s test_data/sample_genome.fasta \
-o tmp/test_output4 -c $NCPU_TO_USE


echo "Running tests 5, detection of LTRs using python wrapper"
./detect_putative_ltr_wrapper.py -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output5 -c $NCPU_TO_USE \
-S 10000000
cat tmp/test_output5_statistics.csv

echo "Running tests 6, filtering gff"
./clean_ltr.R -g tmp/test_output5.gff3 -s test_data/sample_genome.fasta \
-o tmp/test_output6 -c $NCPU_TO_USE

/detect_putative_ltr_wrapper.py -s test_data/sample_genome_part.fasta \
-g test_data/sample_DANTE_part.gff3 -o tmp/test_output7 -c $NCPU_TO_USE \
-S 10000000 -M 2