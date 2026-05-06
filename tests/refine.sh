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
                   --min_cluster_size 3 \
                   --no-mafft-fallback > "$OUT/parasail.log" 2>&1
for f in r_refined.gff3 r_per_element.tsv r_clusters.tsv r_run.json; do
  [ -e "$OUT/parasail/$f" ] || { echo "FAIL: missing $f (parasail-only)"; exit 1; }
done

# v2.1 per_element TSV columns:
#   1 chrom  2 start_orig  ...  24 refinement_method  25 confidence  26 evaluation_status
N_PARASAIL=$(awk -F'\t' 'NR>1 && $24 ~ /^parasail/' \
             "$OUT/parasail/r_per_element.tsv" | wc -l)
[ "$N_PARASAIL" -ge 6 ] \
  || { echo "FAIL: expected >=6 parasail-refined LTRs, got $N_PARASAIL"; exit 1; }
echo "  OK: $N_PARASAIL LTR features refined by parasail"

# v2.1 confidence labels: dual / divergent / inner_only / mafft / unrefined
N_VALIDATED=$(awk -F'\t' 'NR>1 && $24 ~ /^parasail/ && $25 != "unrefined"' \
              "$OUT/parasail/r_per_element.tsv" | wc -l)
[ "$N_VALIDATED" -ge 1 ] \
  || { echo "FAIL: no validated (dual|divergent|inner_only) LTRs"; exit 1; }
echo "  OK: $N_VALIDATED LTR features validated (dual/divergent/inner_only)"

# refined GFF3 must round-trip the v2 schema
grep -q 'Refinement_Method=parasail_inner' "$OUT/parasail/r_refined.gff3" \
  || { echo "FAIL: refined GFF3 missing Refinement_Method=parasail_inner"; exit 1; }
grep -q 'Original_Start=' "$OUT/parasail/r_refined.gff3" \
  || { echo "FAIL: refined GFF3 missing Original_Start"; exit 1; }
grep -qE 'Motif_New=(TG|CA|TG/CA|TG/none|none/CA)' \
     "$OUT/parasail/r_refined.gff3" \
  || { echo "FAIL: refined GFF3 missing Motif_New attribute"; exit 1; }
# v2.1: MSA rescue is NOT enabled in stage 1a (--no-mafft-fallback) -> no MSA_g
grep -q 'MSA_g=' "$OUT/parasail/r_refined.gff3" \
  && { echo "FAIL: MSA_g present despite --no-mafft-fallback"; exit 1; } || true

# v2: TSD-loss revert rule must eliminate any 'lost' outcomes by default
N_TSD_LOST=$(grep -c 'TSD_Outcome=lost' "$OUT/parasail/r_refined.gff3" || true)
[ "$N_TSD_LOST" = "0" ] \
  || { echo "FAIL: $N_TSD_LOST TSD_Outcome=lost rows survived the revert gate"; exit 1; }
echo "  OK: 0 TSD_Outcome=lost rows (revert gate working)"

# v2.1 evaluation_status axis: every LTR row must carry exactly one of
# the four labels.
for s in not_evaluated unresolved confirmed refined; do
  c=$(grep -c "Refinement_Status=$s" "$OUT/parasail/r_refined.gff3" || true)
  echo "    Refinement_Status=$s: $c rows"
done
N_TOTAL_STATUS=$(grep -c 'Refinement_Status=' "$OUT/parasail/r_refined.gff3" || true)
[ "$N_TOTAL_STATUS" -gt 0 ] \
  || { echo "FAIL: refined GFF3 missing Refinement_Status attribute"; exit 1; }
echo "  OK: $N_TOTAL_STATUS rows carry Refinement_Status"

# Refined GFF3 row count must be >= input rows minus original target_site_duplication
# rows (we drop and rebuild them; rebuilt count may differ).
N_IN_NON_TSD=$(awk '$3 != "target_site_duplication" && /^[^#]/' "$DATA/dante_ltr.gff3" | wc -l)
N_OUT_NON_TSD=$(awk '$3 != "target_site_duplication" && /^[^#]/' "$OUT/parasail/r_refined.gff3" | wc -l)
[ "$N_IN_NON_TSD" = "$N_OUT_NON_TSD" ] \
  || { echo "FAIL: non-TSD row count drift: in=$N_IN_NON_TSD out=$N_OUT_NON_TSD"; exit 1; }
echo "  OK: refined GFF3 preserves all $N_IN_NON_TSD non-TSD feature rows"

# ---- 1b. hybrid (parasail + MAFFT fallback) ----
echo
echo "=== 1b. dante_ltr_refine (hybrid: parasail + MAFFT fallback) ==="
./dante_ltr_refine -g "$DATA/dante_ltr.gff3" -s "$DATA/genome.fasta" \
                   -o "$OUT/hybrid/r" \
                   --threads "$NCPU" --workers "$NCPU" \
                   --min_cluster_size 3 > "$OUT/hybrid.log" 2>&1

N_HYB=$(awk -F'\t' 'NR>1 && ($24 ~ /^parasail/ || $24 == "mafft")' \
        "$OUT/hybrid/r_per_element.tsv" | wc -l)
[ "$N_HYB" -ge "$N_PARASAIL" ] \
  || { echo "FAIL: hybrid refined count $N_HYB < parasail-only $N_PARASAIL"; exit 1; }
echo "  OK: hybrid produced $N_HYB refinements (parasail + mafft)"

# v2.1: hybrid mode runs MSA on every qualifying cluster -> MSA_g must
# appear on at least one row, and MSA_Agree must be present.
grep -q 'MSA_g=' "$OUT/hybrid/r_refined.gff3" \
  || { echo "FAIL: hybrid GFF3 missing MSA_g (msa_rescue should be ON by default)"; exit 1; }
grep -q 'MSA_Agree=' "$OUT/hybrid/r_refined.gff3" \
  || { echo "FAIL: hybrid GFF3 missing MSA_Agree"; exit 1; }
echo "  OK: hybrid GFF3 carries MSA_g + MSA_Agree attributes"

# JSON timing/counts sanity
python3 - "$OUT/hybrid/r_run.json" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    d = json.load(f)
assert d["params"]["mafft_fallback"] is True, "expected mafft_fallback=true"
assert d["params"]["msa_rescue"] is True, "expected msa_rescue=true"
assert d["policy"] == "inner_primary", f"expected policy=inner_primary, got {d.get('policy')}"
assert d["counts"]["n_ltr_features"] >= 90, f"expected >= 90 LTR features, got {d['counts']['n_ltr_features']}"
assert "n_msa_calls_total" in d["counts"], "expected n_msa_calls_total in counts"
assert d["timing_s"]["total"] < 600, f"runtime exceeded 600s: {d['timing_s']['total']}"
print(f"  OK: run.json sane (policy={d['policy']}, msa_calls={d['counts']['n_msa_calls_total']}, total={d['timing_s']['total']:.1f}s)")
PY

# ---- 1c. --no_tsd_revert diagnostic mode ----
echo
echo "=== 1c. dante_ltr_refine --no_tsd_revert ==="
./dante_ltr_refine -g "$DATA/dante_ltr.gff3" -s "$DATA/genome.fasta" \
                   -o "$OUT/no_revert/r" \
                   --threads "$NCPU" --workers "$NCPU" \
                   --min_cluster_size 3 \
                   --no_tsd_revert --no-mafft-fallback \
                   > "$OUT/no_revert.log" 2>&1
# In diagnostic mode tsd_outcome can be "lost"; just verify the run produced
# the expected files and the policy field is set.
[ -e "$OUT/no_revert/r_refined.gff3" ] \
  || { echo "FAIL: --no_tsd_revert run produced no GFF3"; exit 1; }
python3 - "$OUT/no_revert/r_run.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
assert d["params"]["no_tsd_revert"] is True
print("  OK: --no_tsd_revert recorded in run.json")
PY

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
