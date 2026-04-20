#!/bin/bash
# tests/fallback.sh — exercises --fallback_mode end-to-end.
# Uses tests/data/fallback_medium (Drapa subset with no aRH; distantly
# related; triggers coarse2/coarse3 recommendation).  Also runs the
# aRH-present branch against the existing sample_DANTE_part data so
# both arms of the data-driven aRH logic are covered.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DRAPA_FA="$ROOT/tests/data/fallback_medium/genome.fasta"
DRAPA_GFF="$ROOT/tests/data/fallback_medium/dante.gff3"
SAMPLE_FA="$ROOT/test_data/sample_genome_part.fasta"
SAMPLE_GFF="$ROOT/test_data/sample_DANTE_part.gff3"
OUT="$ROOT/tmp/tests/fallback"
NCPU="${NCPU:-2}"

for f in "$DRAPA_FA" "$DRAPA_GFF" "$SAMPLE_FA" "$SAMPLE_GFF"; do
  [ -f "$f" ] || { echo "FAIL: $f not found"; exit 1; }
done

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

# ---- 1. Negative control: non-fallback on fallback_medium -> ~0 TEs ----
echo "=== 1. negative control: --fallback_mode=none on fallback_medium ==="
./dante_ltr -g "$DRAPA_GFF" -s "$DRAPA_FA" -o "$OUT/none" -c "$NCPU" \
  > "$OUT/none.log" 2>&1 || true
N_NONE=$(awk -F'\t' '$3=="transposable_element"' "$OUT/none.gff3" 2>/dev/null | wc -l)
echo "  TEs detected (no fallback): $N_NONE"

# Must at least have emitted the recommendation line
grep -q 'Consider re-running with --fallback_mode=coarse2' "$OUT/none.log" \
  || { echo "FAIL: recommendation for coarse2 not printed"; exit 1; }
echo "  OK: recommendation for coarse2 printed"

# ---- 2. coarse2 on fallback_medium -> strictly more TEs than none ----
echo
echo "=== 2. --fallback_mode=coarse2 on fallback_medium ==="
./dante_ltr -g "$DRAPA_GFF" -s "$DRAPA_FA" -o "$OUT/c2" \
            --fallback_mode coarse2 -M 1 -c "$NCPU" > "$OUT/c2.log" 2>&1
N_C2=$(awk -F'\t' '$3=="transposable_element"' "$OUT/c2.gff3" | wc -l)
echo "  TEs detected: $N_C2"
[ "$N_C2" -gt "$N_NONE" ] \
  || { echo "FAIL: coarse2 ($N_C2) not greater than baseline ($N_NONE)"; exit 1; }

# Every classification must be at depth 3 (copia / gypsy)
BAD=$(awk -F'\t' '$3=="transposable_element"' "$OUT/c2.gff3" \
      | grep -oP 'Final_Classification=[^;]+' \
      | awk -F'|' '$0!~/Ty1\/copia$/ && $0!~/Ty3\/gypsy$/ {c++} END{print c+0}')
[ "$BAD" = "0" ] \
  || { echo "FAIL: $BAD TEs at wrong depth in coarse2"; exit 1; }
echo "  OK: every TE classified at depth 3 (copia or gypsy)"

# At least one feature must carry Original_Classification
grep -q 'Original_Classification' "$OUT/c2.gff3" \
  || { echo "FAIL: no Original_Classification in coarse2 output"; exit 1; }
echo "  OK: Original_Classification present on demoted features"

# ---- 3. coarse3 on fallback_medium -> similar or more + chromo rows ----
echo
echo "=== 3. --fallback_mode=coarse3 on fallback_medium ==="
./dante_ltr -g "$DRAPA_GFF" -s "$DRAPA_FA" -o "$OUT/c3" \
            --fallback_mode coarse3 -M 1 -c "$NCPU" > "$OUT/c3.log" 2>&1
N_C3=$(awk -F'\t' '$3=="transposable_element"' "$OUT/c3.gff3" | wc -l)
echo "  TEs detected: $N_C3"
[ "$N_C3" -gt "$N_NONE" ] \
  || { echo "FAIL: coarse3 ($N_C3) not greater than baseline ($N_NONE)"; exit 1; }

# coarse3 depths must be <= 4
BAD=$(awk -F'\t' '$3=="transposable_element"' "$OUT/c3.gff3" \
      | grep -oP 'Final_Classification=[^;]+' \
      | awk -F'|' 'NF > 4 {c++} END {print c+0}')
[ "$BAD" = "0" ] \
  || { echo "FAIL: $BAD TEs exceed depth 4 in coarse3"; exit 1; }
echo "  OK: every TE classified at depth <= 4"

# aRH branch (absent) should have fired with 0%
grep -qE 'aRH check: .* 0\/[0-9]+' "$OUT/c3.log" \
  || grep -qE 'aRH check: .* \(0%\)' "$OUT/c3.log" \
  || { echo "FAIL: aRH-absent path not observed in fallback_medium coarse3"; exit 1; }
echo "  OK: aRH-absent branch exercised (no modified constraints file)"

# ---- 4. coarse3 on sample_DANTE_part -> exercises aRH-PRESENT path ----
echo
echo "=== 4. --fallback_mode=coarse3 on sample_DANTE_part (aRH present) ==="
./dante_ltr -g "$SAMPLE_GFF" -s "$SAMPLE_FA" -o "$OUT/sample_c3" \
            --fallback_mode coarse3 -M 1 -c "$NCPU" > "$OUT/sample_c3.log" 2>&1
N_SC3=$(awk -F'\t' '$3=="transposable_element"' "$OUT/sample_c3.gff3" | wc -l)
echo "  TEs detected: $N_SC3"
[ "$N_SC3" -ge 1 ] \
  || { echo "FAIL: no TEs detected on sample coarse3"; exit 1; }

grep -q 'inserting aRH between RH and INT' "$OUT/sample_c3.log" \
  || { echo "FAIL: aRH insertion did not fire on sample data"; exit 1; }
echo "  OK: aRH inclusion fired and rewrote the non-chromo order"

# The rewritten CSV must exist
find "$(dirname "$OUT/sample_c3")/_fallback_inputs" \
     -name 'constraints_coarse3_arh.csv' 2>/dev/null | head -1 | grep -q . \
  || { echo "FAIL: rewritten constraints file not found"; exit 1; }
echo "  OK: rewritten constraints_coarse3_arh.csv present"

echo
echo "fallback PASSED"
echo "summary:"
echo "  none on fallback_medium : $N_NONE TEs"
echo "  coarse2 on fallback_medium : $N_C2 TEs"
echo "  coarse3 on fallback_medium : $N_C3 TEs"
echo "  coarse3 on sample_part     : $N_SC3 TEs"
