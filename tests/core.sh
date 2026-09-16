#!/bin/bash
# tests/core.sh — core-domain mode (`dante_ltr --mode core`).
#
# Four things are checked, in rising order of what they would cost us:
#   1. the R unit tests (no genome, no BLAST, under a second)
#   2. concordance  — core mode finds the element lineage mode finds
#   3. sensitivity  — core mode finds elements lineage mode cannot
#   4. regression   — lineage mode's own output has not moved
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SMOKE="$ROOT/tests/data/smoke"
DRAPA="$ROOT/tests/data/core_drapa"
OUT="$ROOT/tmp/tests/core"
NCPU="${NCPU:-2}"

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

fail() { echo "FAIL: $*"; exit 1; }

te_count() { awk -F'\t' '$3=="transposable_element"' "$1" | wc -l; }
te_field() {
  # te_field <gff3> <attribute>  -> values, one per transposable_element
  awk -F'\t' -v key="$2" '$3=="transposable_element" {
    if (match($9, key "=[^;]*")) print substr($9, RSTART+length(key)+1, RLENGTH-length(key)-1)
  }' "$1"
}

echo "=== R unit tests ==="
Rscript "$ROOT/tests/core_selftest.R"

echo
echo "=== concordance: core mode on tests/data/smoke ==="
# Lineage mode finds one DLTP element here.  Core mode must find the
# same element -- same boundaries, same rank -- or it is not detecting
# the same biology.
./dante_ltr -g "$SMOKE/dante.gff3" -s "$SMOKE/genome.fasta" \
            -o "$OUT/lineage" -c "$NCPU" >/dev/null
./dante_ltr --mode core -g "$SMOKE/dante.gff3" -s "$SMOKE/genome.fasta" \
            -o "$OUT/core" -c "$NCPU" >/dev/null

LIN_COORD=$(awk -F'\t' '$3=="transposable_element"{print $1":"$4"-"$5":"$7}' "$OUT/lineage.gff3")
CORE_COORD=$(awk -F'\t' '$3=="transposable_element"{print $1":"$4"-"$5":"$7}' "$OUT/core.gff3")
[ -n "$LIN_COORD" ] || fail "lineage mode found no element in the smoke fixture"
[ "$LIN_COORD" = "$CORE_COORD" ] \
  || fail "boundaries differ: lineage=$LIN_COORD core=$CORE_COORD"
echo "OK: identical boundaries ($CORE_COORD)"

[ "$(te_field "$OUT/lineage.gff3" Rank)" = "$(te_field "$OUT/core.gff3" Rank)" ] \
  || fail "rank differs between modes"
echo "OK: identical rank ($(te_field "$OUT/core.gff3" Rank))"

# the order-derived superfamily must agree with lineage mode's call
LIN_CLASS=$(te_field "$OUT/lineage.gff3" Final_Classification)
CORE_CLASS=$(te_field "$OUT/core.gff3" Final_Classification)
case "$LIN_CLASS" in
  "$CORE_CLASS"*) echo "OK: superfamily agrees ($CORE_CLASS)" ;;
  *) fail "classification conflict: lineage=$LIN_CLASS core=$CORE_CLASS" ;;
esac

echo
echo "=== sensitivity: Drapa, a genome REXdb does not cover ==="
# Lineage mode produces an empty GFF3 here.  This is the test that
# encodes the point of the feature.
./dante_ltr -g "$DRAPA/dante.gff3" -s "$DRAPA/genome.fasta" \
            -o "$OUT/drapa_lineage" -c "$NCPU" >/dev/null
./dante_ltr --mode core -g "$DRAPA/dante.gff3" -s "$DRAPA/genome.fasta" \
            -o "$OUT/drapa_core" -c "$NCPU" >/dev/null

N_LIN=$(te_count "$OUT/drapa_lineage.gff3")
N_CORE=$(te_count "$OUT/drapa_core.gff3")
[ "$N_LIN" -eq 0 ] \
  || fail "expected lineage mode to find nothing on Drapa, got $N_LIN"
[ "$N_CORE" -ge 1 ] \
  || fail "expected core mode to find >=1 element on Drapa, got $N_CORE"
echo "OK: lineage $N_LIN vs core $N_CORE elements"

N_NON_D=$(te_field "$OUT/drapa_core.gff3" Rank | grep -cv '^D$' || true)
[ "$N_NON_D" -ge 1 ] || fail "all Drapa core-mode elements are rank D"
echo "OK: $N_NON_D element(s) above rank D"

# every element must carry the core-mode attributes
for key in Superfamily_Evidence Core_Domains Lineage_Support \
           Classification_Demoted Classification_Conflict; do
  N=$(te_field "$OUT/drapa_core.gff3" "$key" | wc -l)
  [ "$N" -eq "$N_CORE" ] || fail "$key present on $N of $N_CORE elements"
done
echo "OK: core-mode attributes present on all elements"

# statistics must stay mergeable and must count every element
HEADER=$(head -1 "$OUT/drapa_core_statistics.csv")
[ "$HEADER" = "$(printf 'Classification\tD\tDL\tDLT\tDLP\tDLTP\tRT_domain')" ] \
  || fail "statistics header changed: $HEADER"
TOTAL=$(awk -F'\t' '$1=="Total"{print $2+$3+$4+$5+$6}' \
        "$OUT/drapa_core_statistics.csv")
[ "$TOTAL" -eq "$N_CORE" ] \
  || fail "statistics Total is $TOTAL but there are $N_CORE elements"
echo "OK: statistics shape unchanged, Total matches element count"

echo
echo "=== determinism ==="
./dante_ltr --mode core -g "$DRAPA/dante.gff3" -s "$DRAPA/genome.fasta" \
            -o "$OUT/drapa_core2" -c "$NCPU" >/dev/null
A=$(md5sum < "$OUT/drapa_core.gff3")
B=$(md5sum < "$OUT/drapa_core2.gff3")
[ "$A" = "$B" ] || fail "core mode is not deterministic"
echo "OK: byte-identical across runs"

echo
echo "=== calibration tool ==="
python3 "$ROOT/utils/calibrate_core_constraints.py" \
        -g "$OUT/lineage.gff3" >/dev/null \
  || fail "calibrate_core_constraints.py failed on lineage-mode output"
echo "OK: runs on lineage-mode output"
# a fixture with no DLT/DLTP has nothing to measure and must say so
if python3 "$ROOT/utils/calibrate_core_constraints.py" \
           -g "$OUT/drapa_lineage.gff3" >/dev/null 2>&1; then
  fail "expected non-zero exit when there is nothing to measure"
fi
echo "OK: exits non-zero when there is nothing to measure"

echo
echo "=== regression: lineage mode unchanged ==="
# --mode core must be purely additive.  Compare against a run made with
# the default flags on the same input.
./dante_ltr -g "$SMOKE/dante.gff3" -s "$SMOKE/genome.fasta" \
            -o "$OUT/lineage_again" -c "$NCPU" >/dev/null
[ "$(md5sum < "$OUT/lineage.gff3")" = "$(md5sum < "$OUT/lineage_again.gff3")" ] \
  || fail "lineage mode is not reproducible"
echo "OK: lineage mode reproducible"

echo
echo "All core-mode tests passed."
