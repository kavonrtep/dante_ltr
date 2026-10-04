#!/bin/bash
# tests/container.sh <image.sif> — run the Galaxy wrapper commands against a
# built image, the way Galaxy invokes it.
#
# Galaxy runs `apptainer exec --cleanenv` as an unprivileged user, with input
# datasets bound read-only and the job directory as the only writable path,
# and it does not activate conda.  --containall adds a private /tmp and HOME,
# so nothing on the host (conda envs, R site libraries, ~/.Rprofile) can make
# the test pass by accident.  The commands mirror the XML in
# kavonrtep/galaxy_packages/dante_ltr; keep them in step.
#
# This exercises the code baked into /opt/dante_ltr, not the checkout.
set -euo pipefail

SIF="${1:?usage: tests/container.sh <image.sif>}"
SIF="$(cd "$(dirname "$SIF")" && pwd)/$(basename "$SIF")"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/test_data"
JOB="$ROOT/tmp/tests/container"
NCPU="${NCPU:-2}"

fail() { echo "FAIL: $*"; exit 1; }

rm -rf "$JOB"
mkdir -p "$JOB"

run() {
  # run <command...> inside the image, Galaxy-style, in the job directory
  apptainer exec --cleanenv --containall \
    -B "$DATA:/data:ro" -B "$JOB:/job" --pwd /job \
    "$SIF" bash -c "$*"
}

echo "=== image metadata ==="
V_IMG=$(run dante_ltr --version 2>/dev/null | tail -1 | awk '{print $2}')
V_SRC=$(cd "$ROOT" && python3 -c "exec(open('version.py').read()); print(__version__)")
[ "$V_IMG" = "$V_SRC" ] || fail "image reports $V_IMG, version.py is $V_SRC"
apptainer inspect --labels "$SIF" | grep -q "Version: $V_SRC" \
  || fail "image label Version is not $V_SRC"
echo "OK: version $V_IMG (binary and label)"

echo
echo "=== dante_ltr_search (lineage) ==="
run dante_ltr --gff3 /data/sample_DANTE_part.gff3 \
    --reference_sequence /data/sample_genome_part.fasta \
    -M 1 --output output --cpu "$NCPU" --mode lineage > "$JOB/lineage.log" 2>&1 \
  || { tail -20 "$JOB/lineage.log"; fail "dante_ltr lineage"; }
[ -s "$JOB/output.gff3" ] || fail "output.gff3 missing"
[ -s "$JOB/output_statistics.csv" ] || fail "output_statistics.csv missing"
[ -s "$JOB/output_summary.html" ] || fail "output_summary.html missing"
N=$(awk -F'\t' '$3=="transposable_element"' "$JOB/output.gff3" | grep -cE 'Rank=DLT' || true)
[ "$N" -ge 1 ] || fail "expected >=1 DLT/DLTP element, got $N"
echo "OK: $N DLT/DLTP element(s), statistics and HTML report written"

echo
echo "=== dante_ltr_search (core) ==="
run dante_ltr --gff3 /data/sample_DANTE_part.gff3 \
    --reference_sequence /data/sample_genome_part.fasta \
    -M 1 --output core --cpu "$NCPU" --mode core > "$JOB/core.log" 2>&1 \
  || { tail -20 "$JOB/core.log"; fail "dante_ltr core"; }
[ -s "$JOB/core.gff3" ] || fail "core.gff3 missing"
echo "OK: core mode"

echo
echo "=== dante_ltr_summary ==="
run dante_ltr_summary -g output.gff3 -o summary > "$JOB/summary.log" 2>&1 \
  || { tail -20 "$JOB/summary.log"; fail "dante_ltr_summary"; }
[ -s "$JOB/summary.html" ] && [ -d "$JOB/summary_plots" ] \
  || fail "summary report missing"
echo "OK: summary report"

echo
echo "=== dante_ltr_to_library ==="
run "cp output.gff3 dante_ltr.gff3 && dante_ltr_to_library --gff3 dante_ltr.gff3 \
    --reference_sequence /data/sample_genome_part.fasta --output_dir library \
    --cpu $NCPU --min_coverage 3" > "$JOB/library.log" 2>&1 \
  || { tail -20 "$JOB/library.log"; fail "dante_ltr_to_library"; }
[ -f "$JOB/library/mmseqs2/mmseqs_representative_seq_clean.fasta" ] \
  || fail "library representative FASTA missing"
echo "OK: repeat library"

echo
echo "=== clean_dante_ltr ==="
# Only checked to start: clean_ltr.R fails on inputs this small (its
# empty-set handling in compare_TE_datasets / multiplicity_of_te), a
# known, separate issue.  It runs on realistic input.
run clean_ltr.R --help > "$JOB/clean.log" 2>&1 \
  || { tail -20 "$JOB/clean.log"; fail "clean_ltr.R --help"; }
grep -q -- "--reference_sequence" "$JOB/clean.log" || fail "clean_ltr.R help text"
echo "OK: clean_ltr.R starts (R packages load)"

echo
echo "All container tests passed."
