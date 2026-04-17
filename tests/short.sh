#!/bin/bash
# tests/short.sh — 1-3 min pipeline check with real positive output.
# Uses the existing test_data/sample_genome_part.fasta (7.5 MB, 3 contigs)
# already in the repo. Asserts non-trivial feature counts and that solo
# LTR detection produces results.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FASTA="$ROOT/test_data/sample_genome_part.fasta"
DANTE="$ROOT/test_data/sample_DANTE_part.gff3"
OUT="$ROOT/tmp/tests/short"
NCPU="${NCPU:-2}"

[ -f "$FASTA" ] || { echo "FAIL: $FASTA not found"; exit 1; }
[ -f "$DANTE" ] || { echo "FAIL: $DANTE not found"; exit 1; }

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

echo "=== dante_ltr on sample_genome_part.fasta ==="
./dante_ltr -g "$DANTE" -s "$FASTA" -o "$OUT/ltr" -c "$NCPU"

# At least one complete element is required for the library build
N_DLTP=$(awk -F'\t' '$3=="transposable_element"' "$OUT/ltr.gff3" \
         | grep -cE 'Rank=(DLTP|DLT)' || true)
[ "$N_DLTP" -ge 1 ] \
  || { echo "FAIL: expected >=1 DLTP or DLT, got $N_DLTP"; exit 1; }
echo "OK: $N_DLTP DLTP/DLT element(s)"

echo
echo "=== dante_ltr_to_library ==="
./dante_ltr_to_library -g "$OUT/ltr.gff3" -s "$FASTA" \
                       -o "$OUT/library" -c "$NCPU"
# library directory has TE_all.fasta
[ -s "$OUT/library/TE_all.fasta" ] \
  || { echo "FAIL: library/TE_all.fasta empty"; exit 1; }
echo "OK: repeat library built"

echo
echo "=== dante_ltr_solo ==="
./dante_ltr_solo -g "$OUT/ltr.gff3" -s "$FASTA" -o "$OUT/solo" -c "$NCPU"

for f in solo_ltr.gff3 solo_ltr_raw.gff3 solo_ltr_statistics.csv; do
  [ -e "$OUT/solo/$f" ] || { echo "FAIL: missing $f"; exit 1; }
done

# Count representatives
N_REP=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr.gff3" | wc -l)
N_RAW=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr_raw.gff3" | wc -l)
echo "OK: $N_RAW raw solo_LTRs collapsed to $N_REP representatives"
[ "$N_REP" -ge 1 ] \
  || { echo "FAIL: expected >=1 solo representative, got $N_REP"; exit 1; }
[ "$N_RAW" -ge "$N_REP" ] \
  || { echo "FAIL: raw ($N_RAW) < representative ($N_REP) is impossible"; exit 1; }

# Every representative must carry a LibraryID
N_MISSING_LIB=$(awk -F'\t' '$3=="solo_LTR" && $9 !~ /LibraryID=LTR_/' \
                "$OUT/solo/solo_ltr.gff3" | wc -l)
[ "$N_MISSING_LIB" = "0" ] \
  || { echo "FAIL: $N_MISSING_LIB solo_LTR(s) missing LibraryID attribute"; exit 1; }
echo "OK: all representatives carry LibraryID"

echo
echo "short PASSED"
