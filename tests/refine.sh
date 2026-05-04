#!/bin/bash
# tests/refine.sh — exercises the per-element boundary refinement
# pipeline added in the parasail-based refinement plan.
#
# Stages:
#   1. dante_ltr_refine on a 1 Mb DANTE_LTR slice (78 TEs across 4
#      lineages) -- both parasail-only and hybrid (parasail + MAFFT
#      fallback) modes.
#   2. utils/build_ltr_library.R --refined_gff3 -- verifies the
#      refined-mode path produces a library and the augmented map TSV.
#   3. dante_ltr_solo --refined_gff3 -- verifies the end-to-end path
#      (refine -> library -> solo detection) produces all expected
#      output files.
#
# Runtime: ~30-60 s on 2 CPUs.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/refine"
OUT="$ROOT/tmp/tests/refine"
NCPU="${NCPU:-2}"

[ -f "$DATA/dante_ltr.gff3" ] || { echo "FAIL: $DATA/dante_ltr.gff3 not found"; exit 1; }
[ -f "$DATA/genome.fasta"   ] || { echo "FAIL: $DATA/genome.fasta not found"; exit 1; }

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$ROOT"

# ---- 1a. parasail-only ----
echo "=== 1a. dante_ltr_refine --no-mafft-fallback ==="
./dante_ltr_refine -g "$DATA/dante_ltr.gff3" -s "$DATA/genome.fasta" \
                   -o "$OUT/parasail/r" \
                   --threads "$NCPU" --workers "$NCPU" \
                   --no-mafft-fallback > "$OUT/parasail.log" 2>&1
for f in r_refined.gff3 r_per_element.tsv r_clusters.tsv r_run.json; do
  [ -e "$OUT/parasail/$f" ] || { echo "FAIL: missing $f (parasail-only)"; exit 1; }
done

N_PARASAIL=$(awk -F'\t' 'NR>1 && $18=="parasail"' "$OUT/parasail/r_per_element.tsv" | wc -l)
[ "$N_PARASAIL" -ge 6 ] \
  || { echo "FAIL: expected >=6 parasail-refined LTRs, got $N_PARASAIL"; exit 1; }
echo "  OK: $N_PARASAIL LTR features refined by parasail"

N_VALIDATED=$(awk -F'\t' 'NR>1 && $18=="parasail" && $19!="low" && $19!="unrefined"' \
              "$OUT/parasail/r_per_element.tsv" | wc -l)
[ "$N_VALIDATED" -ge 1 ] \
  || { echo "FAIL: no validated (high|medium) LTRs"; exit 1; }
echo "  OK: $N_VALIDATED LTR features validated (high or medium confidence)"

# refined GFF3 must round-trip the schema
grep -q 'Refinement_Method=parasail' "$OUT/parasail/r_refined.gff3" \
  || { echo "FAIL: refined GFF3 missing Refinement_Method=parasail"; exit 1; }
grep -q 'Original_Start=' "$OUT/parasail/r_refined.gff3" \
  || { echo "FAIL: refined GFF3 missing Original_Start"; exit 1; }
grep -q 'TG_OK=' "$OUT/parasail/r_refined.gff3" \
  || { echo "FAIL: refined GFF3 missing TG_OK (motif=TG/CA default)"; exit 1; }

# Refined GFF3 row count must match input row count -- pass-through unchanged
N_IN=$(grep -cE '^[^#]' "$DATA/dante_ltr.gff3")
N_OUT=$(grep -cE '^[^#]' "$OUT/parasail/r_refined.gff3")
[ "$N_IN" = "$N_OUT" ] \
  || { echo "FAIL: refined GFF3 row count $N_OUT != input $N_IN"; exit 1; }
echo "  OK: refined GFF3 preserves all $N_IN feature rows"

# ---- 1b. hybrid (parasail + MAFFT fallback) ----
echo
echo "=== 1b. dante_ltr_refine (hybrid: parasail + MAFFT fallback) ==="
./dante_ltr_refine -g "$DATA/dante_ltr.gff3" -s "$DATA/genome.fasta" \
                   -o "$OUT/hybrid/r" \
                   --threads "$NCPU" --workers "$NCPU" > "$OUT/hybrid.log" 2>&1

N_HYB=$(awk -F'\t' 'NR>1 && $18=="parasail"' "$OUT/hybrid/r_per_element.tsv" | wc -l)
[ "$N_HYB" -ge "$N_PARASAIL" ] \
  || { echo "FAIL: hybrid parasail count $N_HYB < parasail-only $N_PARASAIL"; exit 1; }
echo "  OK: hybrid produced $N_HYB parasail refinements"

# JSON timing/counts sanity
python3 - "$OUT/hybrid/r_run.json" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    d = json.load(f)
assert d["params"]["mafft_fallback"] is True, "expected mafft_fallback=true"
assert d["counts"]["n_ltr_features"] >= 90, f"expected >= 90 LTR features, got {d['counts']['n_ltr_features']}"
assert d["timing_s"]["total"] < 600, f"runtime exceeded 600s: {d['timing_s']['total']}"
print(f"  OK: run.json sane (n_ltr_features={d['counts']['n_ltr_features']}, total={d['timing_s']['total']:.1f}s)")
PY

# ---- 1c. parasail-only with motif=none ----
echo
echo "=== 1c. dante_ltr_refine --boundary_motif=none ==="
./dante_ltr_refine -g "$DATA/dante_ltr.gff3" -s "$DATA/genome.fasta" \
                   -o "$OUT/none/r" \
                   --threads "$NCPU" --workers "$NCPU" \
                   --boundary_motif none --no-mafft-fallback \
                   > "$OUT/none.log" 2>&1

# motif=none -> no TG_OK / CA_OK attributes in the refined GFF3
grep -q 'TG_OK=' "$OUT/none/r_refined.gff3" \
  && { echo "FAIL: TG_OK present despite --boundary_motif=none"; exit 1; } || true
grep -q 'Refinement_Method=parasail' "$OUT/none/r_refined.gff3" \
  || { echo "FAIL: motif=none refined GFF3 missing Refinement_Method=parasail"; exit 1; }
echo "  OK: motif=none refined GFF3 has no TG_OK / CA_OK columns"

# ---- 2. build_ltr_library.R --refined_gff3 ----
echo
echo "=== 2. build_ltr_library.R --refined_gff3 ==="
mkdir -p "$OUT/library"
"$ROOT/utils/build_ltr_library.R" \
    --refined_gff3 "$OUT/parasail/r_refined.gff3" \
    -s "$DATA/genome.fasta" \
    -o "$OUT/library/lib" \
    -t "$NCPU" \
    --min_validated_members 4 > "$OUT/library.log" 2>&1

[ -s "$OUT/library/lib_LTR_library.fasta" ] \
  || { echo "FAIL: library FASTA empty"; exit 1; }
[ -s "$OUT/library/lib_LTR_library_map.tsv" ] \
  || { echo "FAIL: library map TSV empty"; exit 1; }

# Map TSV must carry the augmented columns in refined mode
HEAD=$(head -1 "$OUT/library/lib_LTR_library_map.tsv")
echo "$HEAD" | grep -q 'low_confidence' \
  || { echo "FAIL: refined map TSV missing low_confidence column"; exit 1; }
echo "$HEAD" | grep -q 'consensus_built_from' \
  || { echo "FAIL: refined map TSV missing consensus_built_from column"; exit 1; }
echo "$HEAD" | grep -q 'n_validated' \
  || { echo "FAIL: refined map TSV missing n_validated column"; exit 1; }
echo "  OK: map TSV has refined-mode columns"

N_VALIDATED_CL=$(awk -F'\t' 'NR>1 && $5=="validated"' \
                 "$OUT/library/lib_LTR_library_map.tsv" | wc -l)
echo "  OK: $N_VALIDATED_CL cluster(s) built consensus from validated subset"

# ---- 3. dante_ltr_solo --refined_gff3 (end-to-end) ----
echo
echo "=== 3. dante_ltr_solo --refined_gff3 ==="
./dante_ltr_solo -g "$DATA/dante_ltr.gff3" \
                 --refined_gff3 "$OUT/parasail/r_refined.gff3" \
                 -s "$DATA/genome.fasta" \
                 -o "$OUT/solo" -c "$NCPU" \
                 --min_validated_members 4 > "$OUT/solo.log" 2>&1

for f in solo_ltr.gff3 solo_ltr_raw.gff3 solo_ltr_statistics.csv \
         reference_input_refined.gff3; do
  [ -e "$OUT/solo/$f" ] || { echo "FAIL: missing solo output $f"; exit 1; }
done
[ -s "$OUT/solo/library/ltr_lib_LTR_library_map.tsv" ] \
  || { echo "FAIL: missing solo library map"; exit 1; }
head -1 "$OUT/solo/library/ltr_lib_LTR_library_map.tsv" | grep -q 'low_confidence' \
  || { echo "FAIL: solo library map missing refined-mode columns"; exit 1; }
echo "  OK: dante_ltr_solo --refined_gff3 produced refined library + solo outputs"

echo
echo "refine PASSED"
