#!/bin/bash
# tests/long.sh — full pipeline on the g1 test data (~10-30 min).
# Gate for release tags. Uses loose count ranges, not exact numbers,
# so minor version-to-version drift doesn't fail CI.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Default: the committed sample_genome.fasta + sample_DANTE.gff3.
# Set LONG_FASTA / LONG_DANTE to point at a bigger dataset
# (e.g. test_data/g1.fasta + test_data/g1_dante_ltr.gff3) for deeper
# local testing; those are not in the repo.
FASTA="${LONG_FASTA:-$ROOT/test_data/sample_genome.fasta}"
DANTE="${LONG_DANTE:-$ROOT/test_data/sample_DANTE.gff3}"
OUT="$ROOT/tmp/tests/long"
NCPU="${NCPU:-4}"

[ -f "$FASTA" ] || { echo "FAIL: $FASTA not found"; exit 1; }
[ -f "$DANTE" ] || { echo "FAIL: $DANTE not found"; exit 1; }

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

# Determine whether the DANTE gff is raw DANTE output or already dante_ltr output.
if grep -q $'\ttransposable_element\t' "$DANTE"; then
  echo "=== input already contains dante_ltr features, skipping dante_ltr ==="
  LTR_GFF="$DANTE"
else
  echo "=== dante_ltr ==="
  ./dante_ltr -g "$DANTE" -s "$FASTA" -o "$OUT/ltr" -c "$NCPU"
  LTR_GFF="$OUT/ltr.gff3"
fi

N_TE=$(awk -F'\t' '$3=="transposable_element"' "$LTR_GFF" | wc -l)
echo "OK: $N_TE transposable_element feature(s)"
[ "$N_TE" -ge 50 ] \
  || { echo "FAIL: expected >=50 TEs, got $N_TE"; exit 1; }

echo
echo "=== dante_ltr_solo ==="
./dante_ltr_solo -g "$LTR_GFF" -s "$FASTA" -o "$OUT/solo" -c "$NCPU"

N_REP=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr.gff3" | wc -l)
N_SL=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr.gff3" \
       | grep -c 'Rank=SL;' || true)
echo "OK: $N_REP representatives ($N_SL with confirmed TSD)"
[ "$N_REP" -ge 50 ] \
  || { echo "FAIL: expected >=50 representatives, got $N_REP"; exit 1; }
[ "$N_SL"  -ge  1 ] \
  || { echo "FAIL: expected >=1 SL (TSD-confirmed), got $N_SL"; exit 1; }

echo
echo "long PASSED"
