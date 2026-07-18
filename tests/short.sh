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

# Count representatives (after collapse + TE-fragment partition)
N_REP=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr.gff3" | wc -l)
N_RAW=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr_raw.gff3" | wc -l)
echo "OK: $N_RAW raw solo_LTRs -> $N_REP solo representatives (TE fragments split off separately)"
[ "$N_RAW" -ge "$N_REP" ] \
  || { echo "FAIL: raw ($N_RAW) < representative ($N_REP) is impossible"; exit 1; }

# Every representative must carry a LibraryID
N_MISSING_LIB=$(awk -F'\t' '$3=="solo_LTR" && $9 !~ /LibraryID=LTR_/' \
                "$OUT/solo/solo_ltr.gff3" | wc -l)
[ "$N_MISSING_LIB" = "0" ] \
  || { echo "FAIL: $N_MISSING_LIB solo_LTR(s) missing LibraryID attribute"; exit 1; }
echo "OK: all representatives carry LibraryID"

# Auto-refinement: dante_ltr_solo should have produced a refined GFF3
# under <out>/refined/ (since the input was an unrefined DANTE_LTR GFF3).
[ -e "$OUT/solo/refined/sample_refined.gff3" ] \
  || { echo "FAIL: dante_ltr_solo did not auto-produce refined/sample_refined.gff3"; exit 1; }
echo "OK: dante_ltr_solo auto-refined the input"

# LibraryConfidence propagation onto solo_LTR features
N_LIBCONF=$(awk -F'\t' '$3=="solo_LTR" && $9 ~ /LibraryConfidence=/' \
            "$OUT/solo/solo_ltr.gff3" | wc -l)
[ "$N_LIBCONF" -ge 1 ] \
  || { echo "FAIL: no solo_LTR features carry LibraryConfidence attribute"; exit 1; }
echo "OK: $N_LIBCONF solo_LTR(s) carry LibraryConfidence"

# TE-fragment split: solo_ltr_te_fragments.gff3 must exist; every entry
# must be SL_noTSD with at least one positive junction; the main file
# must be free of such entries.
[ -e "$OUT/solo/solo_ltr_te_fragments.gff3" ] \
  || { echo "FAIL: missing solo_ltr_te_fragments.gff3"; exit 1; }
N_FRAG=$(awk -F'\t' '$3=="solo_LTR"' "$OUT/solo/solo_ltr_te_fragments.gff3" | wc -l)
N_BAD_FRAG=$(awk -F'\t' '$3=="solo_LTR" && ($9 !~ /Rank=SL_noTSD/ || \
              ($9 !~ /UTR5_junction=positive/ && \
               $9 !~ /PPT_junction=positive/  && \
               $9 !~ /PBS_check=positive/))' \
              "$OUT/solo/solo_ltr_te_fragments.gff3" | wc -l)
[ "$N_BAD_FRAG" = "0" ] \
  || { echo "FAIL: $N_BAD_FRAG entry/ies in te_fragments don't match the rule"; exit 1; }
N_LEAK=$(awk -F'\t' '$3=="solo_LTR" && /Rank=SL_noTSD/ && \
              (/UTR5_junction=positive/ || /PPT_junction=positive/ || /PBS_check=positive/)' \
         "$OUT/solo/solo_ltr.gff3" | wc -l)
[ "$N_LEAK" = "0" ] \
  || { echo "FAIL: $N_LEAK probable TE fragment(s) leaked into main solo_ltr.gff3"; exit 1; }
echo "OK: $N_FRAG TE fragment(s) partitioned cleanly out of main solo file"

echo
echo "=== dante_ltr_refine (parasail anchored extension) ==="
./dante_ltr_refine -g "$OUT/ltr.gff3" -s "$FASTA" \
                   -o "$OUT/refined/sample" \
                   --threads "$NCPU" --workers "$NCPU" --no-mafft-fallback

for f in sample_refined.gff3 sample_per_element.tsv sample_clusters.tsv sample_run.json; do
  [ -e "$OUT/refined/$f" ] || { echo "FAIL: missing refine output $f"; exit 1; }
done

# refined GFF3 should have at least one Refinement_Method attribute
N_REFINED_ATTR=$(grep -c 'Refinement_Method=' "$OUT/refined/sample_refined.gff3" || true)
[ "$N_REFINED_ATTR" -ge 1 ] \
  || { echo "FAIL: refined GFF3 has no Refinement_Method attributes"; exit 1; }
echo "OK: refined GFF3 has $N_REFINED_ATTR Refinement_Method attributes"

# per-element TSV header sanity
head -1 "$OUT/refined/sample_per_element.tsv" | grep -q refinement_method \
  || { echo "FAIL: per_element.tsv missing refinement_method header"; exit 1; }

echo
echo "=== dante_ltr_solo chunked path ==="
# Regression guard for the large-genome branch (_run_chunked). It is only
# taken when the reference exceeds --max_chunk_size AND has >1 sequence, so
# the run above (default -S 100 Mbp) never reaches it -- a break there once
# went unnoticed by the whole suite. Forced here with a small -S on the same
# 3-contig fixture. The already-refined GFF3 is reused so the run does not
# pay for refinement a second time.
./dante_ltr_solo -g "$OUT/solo/refined/sample_refined.gff3" -s "$FASTA" \
                 -o "$OUT/solo_chunked" -c "$NCPU" -S 1000000 >/dev/null

grep -q "Running analysis in chunks" "$OUT/solo_chunked/dante_ltr_solo.log" \
  || { echo "FAIL: chunked path was not taken -- the test no longer guards it"; exit 1; }
N_CHUNK_FILES=$(grep -c "Detecting solo LTRs on chunk_" \
                "$OUT/solo_chunked/dante_ltr_solo.log" || true)
[ "$N_CHUNK_FILES" -ge 2 ] \
  || { echo "FAIL: expected >=2 chunks, got $N_CHUNK_FILES"; exit 1; }
echo "OK: chunked path ran over $N_CHUNK_FILES chunks"

# Chunking must not change the calls: coordinates have to survive the split
# and the remap back. Only the feature IDs legitimately differ -- they carry
# a per-chunk serial and the chunk name (ctg16_0) rather than the contig --
# and the remap drops the trailing ';' that the unchunked path leaves in
# place, so both are normalized away before comparing.
norm_gff3() {
  grep -v '^#' "$1" | sed -E 's/(ID|Parent)=[^;]*;?//g; s/;$//' | sort
}
for f in solo_ltr.gff3 solo_ltr_raw.gff3 solo_ltr_te_fragments.gff3; do
  [ -e "$OUT/solo_chunked/$f" ] \
    || { echo "FAIL: chunked run missing $f"; exit 1; }
  diff <(norm_gff3 "$OUT/solo/$f") <(norm_gff3 "$OUT/solo_chunked/$f") >/dev/null \
    || { echo "FAIL: chunked $f differs from unchunked beyond feature IDs"; exit 1; }
done
echo "OK: chunked output matches unchunked on all three GFF3s"

echo
echo "short PASSED"
