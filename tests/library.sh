#!/bin/bash
# tests/library.sh — the repeat-library annotation policy.
#
# Two things are checked, and the first matters more than the second:
#   1. --annotation_conflict strict is byte-identical to the historical code
#   2. the policy itself behaves as specified, driven directly with synthetic
#      clusters
#
# The end-to-end nested path is deliberately NOT asserted on a real genome
# here: no fixture in this repository is large enough to form a cluster whose
# members disagree, and on the genome that motivated the feature both counting
# conventions give the same answer. The unit tests are the real coverage.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/test_data"
OUT="$ROOT/tmp/tests/library"
NCPU="${NCPU:-2}"

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$ROOT"
fail() { echo "FAIL: $*"; exit 1; }

echo "=== policy unit tests ==="
Rscript "$ROOT/tests/library_policy_selftest.R"

echo
echo "=== CLI surface ==="
./dante_ltr_to_library --help | grep -q -- '--annotation_conflict' \
  || fail "--annotation_conflict missing from --help"
./dante_ltr_to_library --help | grep -q -- '--lineage_promotion_min_share' \
  || fail "--lineage_promotion_min_share missing from --help"
./utils/mmseq_clustering.R --help 2>&1 | grep -q -- '--annotation_conflict' \
  || fail "--annotation_conflict missing from mmseq_clustering.R"
echo "OK: options exposed on both entry points"

echo
echo "=== strict is the default, and is byte-identical ==="
./dante_ltr -g "$DATA/sample_DANTE_part.gff3" -s "$DATA/sample_genome_part.fasta" \
            -o "$OUT/ltr" -c "$NCPU" >/dev/null 2>&1
./dante_ltr_to_library -g "$OUT/ltr.gff3" -s "$DATA/sample_genome_part.fasta" \
            -o "$OUT/lib_default" -c "$NCPU" > "$OUT/default.log" 2>&1
./dante_ltr_to_library -g "$OUT/ltr.gff3" -s "$DATA/sample_genome_part.fasta" \
            -o "$OUT/lib_strict" -c "$NCPU" -a strict > "$OUT/strict.log" 2>&1

REP=mmseqs2/mmseqs_representative_seq_clean.fasta
RMC=mmseqs2/mmseqs_representative_seq_clean_rm_compatible.fasta
for f in "$REP" "$RMC"; do
  [ -s "$OUT/lib_default/$f" ] || fail "missing $f"
  a=$(md5sum < "$OUT/lib_default/$f"); b=$(md5sum < "$OUT/lib_strict/$f")
  [ "$a" = "$b" ] || fail "explicit --annotation_conflict strict differs from the default"
done
echo "OK: default == strict, $(grep -c '^>' "$OUT/lib_default/$REP") sequences"

# The summary line has to be present and has to account for every cluster.
grep -q "annotation policy strict:" "$OUT/default.log" \
  || fail "summary line not printed"
python3 - "$OUT/default.log" <<'PY'
import re, sys
m = re.search(r"annotation policy \w+: (\d+) clusters \| kept (\d+) "
              r"\(majority (\d+), lca (\d+), recovered (\d+), promoted (\d+)\) "
              r"\| dropped (\d+)", open(sys.argv[1]).read())
assert m, "summary line did not parse"
tot, kept, maj, lca, rec, pro, drop = map(int, m.groups())
assert maj + lca + rec + pro == kept, "kept breakdown does not sum"
assert kept + drop == tot, "kept + dropped != total"
print("OK: summary line accounts for all %d clusters" % tot)
PY

echo
echo "=== nested runs end-to-end and does not crash ==="
./dante_ltr_to_library -g "$OUT/ltr.gff3" -s "$DATA/sample_genome_part.fasta" \
            -o "$OUT/lib_nested" -c "$NCPU" -a nested > "$OUT/nested.log" 2>&1
[ -s "$OUT/lib_nested/$REP" ] || fail "nested produced no library"
grep -q "annotation policy nested:" "$OUT/nested.log" || fail "nested summary missing"
echo "OK: $(grep -c '^>' "$OUT/lib_nested/$REP") sequences"

# On this fixture every cluster is uniform, so nested must agree with strict.
# Where it would differ there is a disagreeing cluster, which this genome has
# none of -- so a difference here means the policy fired when it should not.
a=$(md5sum < "$OUT/lib_strict/$REP"); b=$(md5sum < "$OUT/lib_nested/$REP")
[ "$a" = "$b" ] || fail "nested changed a library with no disagreeing clusters"
echo "OK: nested is a no-op where no cluster disagrees"

echo
echo "=== rejects a bad policy name ==="
if ./utils/mmseq_clustering.R -f /dev/null -o "$OUT/bad" -a bogus >/dev/null 2>&1; then
  fail "an invalid --annotation_conflict was accepted"
fi
echo "OK: invalid policy rejected"

echo
echo "All library-policy tests passed."
