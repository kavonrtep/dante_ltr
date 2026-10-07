#!/bin/bash
# tests/smoke.sh — fast sanity check (< 30 s).
# Verifies every CLI responds and the dante_ltr pipeline runs on a tiny
# 40 kb genome slice with one annotated DLTP element.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/smoke"
OUT="$ROOT/tmp/tests/smoke"
NCPU="${NCPU:-2}"

rm -rf "$OUT"
mkdir -p "$OUT"

cd "$ROOT"

echo "=== CLI --version / --help ==="
./dante_ltr            --version
./dante_ltr_solo       --help    > /dev/null
./dante_ltr_to_library --help    > /dev/null
./dante_ltr_summary    --help    > /dev/null

echo
echo "=== dante_ltr on 40 kb subset ==="
./dante_ltr -g "$DATA/dante.gff3" -s "$DATA/genome.fasta" \
            -o "$OUT/ltr" -c "$NCPU" >/dev/null

# Assertions
[ -s "$OUT/ltr.gff3" ] || { echo "FAIL: empty ltr.gff3"; exit 1; }
head -1 "$OUT/ltr.gff3" | grep -q '^##gff-version 3' \
  || { echo "FAIL: missing GFF3 header"; exit 1; }
N_TE=$(awk -F'\t' '$3=="transposable_element"' "$OUT/ltr.gff3" | wc -l)
[ "$N_TE" -ge 1 ] \
  || { echo "FAIL: expected >=1 transposable_element, got $N_TE"; exit 1; }
echo "OK: $N_TE transposable_element feature(s)"

echo
echo "=== opposite-strand domain inside an element (issue #15) ==="
# Must complete with -M 1 instead of aborting in get_ranges_right(); the
# output has no complete element, so this also covers the empty summary.
OSO="$ROOT/tests/data/opposite_strand_overlap"
./dante_ltr -g "$OSO/dante.gff3" -s "$OSO/genome.fasta" \
            -o "$OUT/oso" -c "$NCPU" -M 1 >/dev/null 2>&1 \
  || { echo "FAIL: dante_ltr aborted on opposite_strand_overlap"; exit 1; }
[ -s "$OUT/oso_summary.html" ] \
  || { echo "FAIL: no summary report for output without complete elements"; exit 1; }
echo "OK: completes, summary written"

echo
echo "=== dante_ltr_solo pipeline ==="
./dante_ltr_solo -g "$OUT/ltr.gff3" -s "$DATA/genome.fasta" \
                 -o "$OUT/solo" -c "$NCPU" >/dev/null

for f in solo_ltr.gff3 solo_ltr_raw.gff3 solo_ltr_statistics.csv \
         solo_ltr_raw_statistics.csv; do
  [ -e "$OUT/solo/$f" ] || { echo "FAIL: missing $f"; exit 1; }
done
head -1 "$OUT/solo/solo_ltr.gff3" | grep -q '^##gff-version 3' \
  || { echo "FAIL: solo_ltr.gff3 missing GFF3 header"; exit 1; }
echo "OK: solo pipeline produced all expected output files"

echo
echo "smoke PASSED"
