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
./dante_ltr -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output1 -c $NCPU_TO_USE


cat tmp/test_output1_statistics.csv

echo "Running tests 2, create library"
./dante_ltr_to_library -s test_data/sample_genome.fasta \
-g tmp/test_output1.gff3 -o tmp/test_output2 -c $NCPU_TO_USE

echo "Running tests 3, detection of LTRs, allow missing domains"
./dante_ltr -s test_data/sample_genome.fasta \
-g test_data/sample_DANTE.gff3 -o tmp/test_output3 -c $NCPU_TO_USE -M 2


cat tmp/test_output3_statistics.csv

echo "Running tests 4, create library"
./dante_ltr_to_library -s test_data/sample_genome.fasta \
-g tmp/test_output3.gff3 -o tmp/test_output4 -c $NCPU_TO_USE


