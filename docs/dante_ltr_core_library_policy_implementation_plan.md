# Library annotation policy for core mode — implementation plan (v1)

**Request:** `docs/dante_ltr_core_library_policy_request.md` (R1–R5).
**Scope:** `dante_ltr_to_library` / `utils/mmseq_clustering.R` against current
`main` (0.6.0.0). No change to detection, to `--mode core`, or to any other
entry point.
**Testing discipline:** the default path must stay **byte-identical**. That is
the gate for every commit, not just the last one. `smoke` + `short` + a new
`library` target pass locally before each commit; no push until the whole plan
lands.

---

## 1. Files touched

```
utils/library_policy.R          [new]      the keep/label decision, as a pure function
utils/mmseq_clustering.R        [modified] call it; add the three options; R5 log line
dante_ltr_to_library            [modified] plumb the three options through
tests/library_policy_selftest.R [new]      unit tests for the policy, no mmseqs needed
tests/library.sh                [new]      byte-identity + nested end-to-end
tests.sh                        [modified] 'library' target
.github/workflows/tests.yml     [modified] run it
README.md, changelog.md, version.py        [modified] release 0.6.1.0
```

Not touched: `utils/extract_fasta.R`, anything under detection.

---

## 2. The baseline, and the one convention that must be preserved

Resolved by the requester: the earlier 393/397 gap was a tie-break difference,
not element collapse. `which.max(table(x))` breaks ties by taking the **first
label in sorted order**; `Counter.most_common` takes **insertion order**. Three
clusters differed. A faithful port reproduces the shipped 397 exactly.

That makes the tie-break a **requirement on this work, not a footnote**. The
policy function must keep `which.max(table(x))` semantics, because the tie-break
is the difference between a deterministic library and one that depends on record
order — the same class of problem the canonical sort at
`mmseq_clustering.R:90` already exists to prevent. Step 1 must not "tidy" it
into a `which.max` over a differently-ordered table, and Step 5 tests it
directly with a deliberate two-way tie.

Corrected baseline, all four categories over the same pool:

```
                       pool   kept                    dropped
strict (windows)        464    397 = 355 + 42          67 = 62 chain + 5 conflict
nested (elements)       464    459 = 355 + 42 + 62      5 = conflict
delta                         +62 sequences, +68,162 bp  (+16 %)
```

**Element collapse changes no keep/drop decision on this genome** — both
conventions partition the 464 identically. Two consequences for the plan:

- the 459 target is *insensitive* to R3, so it does not validate R3 at all;
- R3 is therefore only reachable by unit test (§4 Step 5), which is now its
  sole coverage rather than a convenience.

## 3. Design: the decision becomes a pure function

Today the policy is four lines of vectorised code fused into the script body
(`mmseq_clustering.R:129-139`), which cannot be tested without running mmseqs
over a real library. Extracting it is what makes R2–R4 verifiable at all.

`utils/library_policy.R` exposes one function:

```r
classify_cluster(labels, elements, proportion_min, policy,
                 promote_min_elements, promote_min_share)
  -> list(keep =, label =, outcome =)     # outcome: "majority" | "lca" |
                                          # "recovered" | "promoted" | "dropped"
```

`labels` and `elements` are parallel vectors, one entry per cluster member
(window). Everything the policy needs is local to one cluster, so the function
is total and trivially testable.

### 3.1 `policy = "strict"` — today's rule, restated

Counts **windows**, exactly as now:

```
prop  <- max(table(labels)) / length(labels)
main  <- names(which.max(table(labels)))
label <- if (prop > proportion_min) main else resolve_name(unique(labels))
keep  <- label == main
```

Implemented by calling the existing `resolve_name()` unchanged, so the strict
branch is the current code moved, not rewritten.

### 3.2 `policy = "nested"` — R2, R3, R4

Counts **distinct elements** (R3). One element contributes one vote regardless
of how many 1 kb windows it was cut into:

```
el     <- unique(data.frame(element, label))      # one row per element
n      <- nrow(el)
prop   <- max(table(el$label)) / n
main   <- names(which.max(table(el$label)))
chain  <- is_single_chain(unique(el$label))
```

`is_single_chain()` (R2): sort the distinct labels by depth; every label must
be a prefix of the next **at `|` boundaries**. Then:

```
keep <- (prop > proportion_min) || chain
base <- if (prop > proportion_min) main else shallowest(chain)
```

and promotion (R4) is evaluated in the **chain branch only** — see Step 4 for
why it cannot fire under a >95 % majority:

```
deep   <- deepest(chain)
n_deep <- number of distinct elements carrying `deep`
label  <- if (n_deep >= promote_min_elements &&
              n_deep / n >= promote_min_share)  deep  else  base
```

Sibling and cross-superfamily mixes fail `is_single_chain()` and are dropped —
with one exception that needs a decision, in §3.5.

### 3.3 Why the chain test rather than "is the LCA present"

R2 asks for the stricter form and it is worth honouring: `is_single_chain()`
can never place two named lineages under one representative, whereas "LCA is
present" would accept `{gypsy, chromovirus|Tekay, chromovirus|Reina}`. On the
Drapa genome the two agree, so this costs nothing and removes a failure mode.

The chain test also lets the nested path avoid `resolve_name()` entirely: the
LCA of a chain *is* its shallowest member, already in hand. That matters
because `resolve_name()` takes `max()` over the positions at which all labels
agree, which would be wrong for a non-prefix-consistent vocabulary. Not a bug
today, but not something to build new behaviour on.

### 3.4 Invariants to assert in code

Both are cheap and are what makes promotion safe (R4):

```r
stopifnot(is_descendant_or_equal(label, base))   # refinement, never reclassification
stopifnot(!promoted || chain)                    # never with siblings present
```

### 3.5 `nested` is not automatically a superset of `strict`

R2 says sibling mixes "stay dropped, exactly as today", but today they are not
always dropped. `strict` keeps any cluster whose LCA equals its majority label,
and that can include a sibling mix:

```
{chromovirus|Tekay, chromovirus|Reina, chromovirus}   majority = chromovirus
  strict : LCA = chromovirus = majority  -> KEPT, labelled chromovirus
  nested : not a single chain            -> DROPPED
```

So a cluster in the library today could vanish under `nested`. On this genome it
does not happen — the requester's 459 = 355 + 42 + 62 means all 42
`LCA == majority` clusters are also chains — but that is a property of the data,
not of the rule.

**Decided: follow R2 literally.** The keep rule is

```
keep <- (prop > proportion_min) || chain
```

with no carry-over of strict's `lca == main` term. `nested` is therefore *not*
guaranteed to be a superset of `strict`: a sibling-plus-parent cluster that
`strict` keeps is dropped. Identical outcome on this genome; the consequence is
a documentation obligation, not a code one —

- the README and the changelog must state that `nested` may **remove** library
  sequences as well as add them, and
- `tests/library_policy_selftest.R` pins the behaviour with the
  sibling-plus-parent case, so the asymmetry is deliberate and visible rather
  than discovered later on someone's genome.

---

## 4. Implementation order

### Step 1 — extract the policy, no behaviour change

Move the decision into `utils/library_policy.R` with `policy = "strict"` only,
and have `mmseq_clustering.R` call it per cluster. Add
`tests/library_policy_selftest.R` covering the strict cases.

This is the risky commit, because it converts vectorised code into a per-cluster
call. **Gate: `short` output MD5-identical.**

Fold in one latent bug here, since this commit is already the byte-identity
checkpoint:

```r
all_names <- sapply(annot_in_clusters, function(x) unique(x))
```

`sapply` simplifies to a **matrix** when every cluster happens to have the same
number of distinct labels, after which `sapply(all_names, resolve_name)` returns
one value per label instead of per cluster and the `ifelse` silently recycles.
Verified in R: two clusters of two labels each give a 2×2 matrix and a
length-4 result where 2 is expected. It cannot fire on a ragged real run, but it
fires on exactly the kind of small fixture Step 5 adds. The per-cluster function
removes it structurally.

**Commit:** `refactor: extract the library annotation policy into a function`

### Step 2 — R1, the flag, still strict-only

`--annotation_conflict {strict,nested}` on both `dante_ltr_to_library` and
`mmseq_clustering.R`, default `strict`, plus
`--lineage_promotion_min_elements` (2) and `--lineage_promotion_min_share`
(0.25). `nested` accepted and parsed but not yet implemented — keeps the CLI
change separable from the behaviour change.

Note while here: `dante_ltr_to_library`'s existing `-p/--proportion_min` has no
`type=`, so it arrives as a string. Harmless (it is interpolated into a shell
command) but the new numeric options should declare `type=float` / `type=int`.

**Gate: `short` MD5-identical; `--help` shows the new options.**

**Commit:** `library: add --annotation_conflict, defaulting to current behaviour`

### Step 3 — R2 + R3, recover ancestor chains

Implement the `nested` branch: element-level counting and `is_single_chain()`.
Promotion not yet wired, so a recovered cluster takes the LCA.

**Gate:** strict still MD5-identical; on the request's genome, nested keeps
355 + 42 + 62 = 459 clusters, dropping only the 5 cross-superfamily ones
(+62 sequences, +68,162 bp against strict's 397).

**Commit:** `library: recover ancestor-chain clusters under --annotation_conflict nested`

### Step 4 — R4, promotion

Add the two thresholds and the two `stopifnot` invariants.

Promotion is evaluated **only in the chain branch**, and deliberately skipped
where the majority rule fired. That is sound rather than incidental: a >95 %
majority leaves every other label under 5 %, which cannot clear a 25 % share, and
a deepest label that *is* the majority is already the cluster's label. Measured
by prior outcome: majority 0, `LCA == majority` 10, recovered 59. Write the
argument into the code comment so a later reader does not "fix" the omission.

**Gate:** 69 clusters promoted on the request's genome — 50 to `chromovirus`,
19 to `Chlamyvir`, split 0 / 10 / 59 by prior outcome; and the promoted label is
a descendant of the LCA in every case.

**Commit:** `library: promote clusters to the deepest supported label`

### Step 5 — R5 and tests

One summary line after clustering:

```
annotation policy nested: 464 clusters | kept 397 | recovered 62 | promoted 69 | dropped 5
```

`tests/library.sh`:
- byte-identity: `dante_ltr_to_library` on `tests/data/smoke` output with
  default flags, MD5 against a recorded baseline;
- `tests/short.sh` already builds a library and has a determinism guard —
  extend it rather than duplicating the fixture;
- a `nested` end-to-end on a **synthetic** cluster fixture (below);
- the policy unit tests from Step 1, extended to R2–R4.

Fixture: the real fixtures are too small to form clusters of ≥3 windows
(`smoke` has one element), so the nested cases are covered by driving
`classify_cluster()` directly with synthetic label/element vectors. That is
also the only way to test the counting rule of R3 — window vs element — without
constructing a genome whose elements differ in length by 8×.

Cases: chain of 2 and of 3; single label; siblings; cross-superfamily; chain
where promotion fires; chain where it fails on `min_elements`; chain where it
fails on `min_share`; the same cluster counted by windows vs elements giving
different majorities; **a deliberate two-way tie**, asserting the sorted-order
tie-break of §2; and a sibling pair whose parent holds the majority (§6).

**Commit:** `library: summary line and tests for the annotation policy`

### Step 6 — docs and release

README section under the core-mode text, changelog, `version.py` → 0.6.1.0.
State plainly that `nested` is for core-mode output and that it is opt-in
because lineage mode can also emit internal-node labels (R1's reasoning).

**Commit:** `docs: document --annotation_conflict; release 0.6.1.0`

---

## 5. Local verification gates

| gate | command | must show |
|---|---|---|
| default unchanged | `./tests.sh short` | library MD5 equal to pre-change |
| policy unit tests | `Rscript tests/library_policy_selftest.R` | all pass |
| new target | `./tests.sh library` | pass |
| nothing else moved | `./tests.sh smoke fallback core refine` | pass |

Capture the baseline **before Step 1**:

```bash
./dante_ltr -g test_data/sample_DANTE_part.gff3 -s test_data/sample_genome_part.fasta \
            -o /tmp/base -c 4
./dante_ltr_to_library -g /tmp/base.gff3 -s test_data/sample_genome_part.fasta \
            -o /tmp/base_lib -c 4
md5sum /tmp/base_lib/mmseqs2/mmseqs_representative_seq_clean*.fasta > /tmp/lib_baseline.md5
```

---

## 6. Risks

| risk | mitigation |
|---|---|
| Step 1's refactor changes strict output | it is the whole point of the gate; if MD5 moves, revert and do the extraction in smaller pieces |
| Element-id parsing wrong | member names are `<seqid>_<start>_<end>#<classification>_sliding:<w0>-<w1>`; confirmed no `#` or `_sliding` occurs inside an id on the Drapa run. Assert it at runtime and fail loudly rather than mis-grouping |
| `min_coverage` unit (see §7) | left on windows, so strict and nested see the same cluster set and only the decision differs |
| Promotion mislabels a cluster | the two `stopifnot` invariants; promotion can only move *down one chain* |
| mmseqs order-sensitivity confounds the comparison | already handled — `mmseq_clustering.R` sorts canonically, and `short.sh` guards it. Do not touch that sort |
| Tie-break silently changed during the refactor | §2: `which.max(table(x))` is sorted-order, not first-seen. A deliberate two-way tie is in the Step 5 unit tests |
| `nested` drops a cluster `strict` keeps | accepted by decision (§3.5) — sibling-plus-parent mixes are dropped. Pinned by a unit test and called out in the release note; does not occur on the request's genome |

---

## 7. Open questions for the requester

**Q1 — `-m/--min_coverage` stays on windows. Answered and quantified.**
17 of the 464 clusters have fewer than 3 distinct elements, and they are in the
library under `strict` today, so this is a pre-existing property rather than
something `nested` introduces. Changing the threshold would move the pool and
invalidate the 459 target, so it is out of scope here and recorded as a separate
decision.

**Q2 — answered: follow R2 literally.** `nested` may drop a sibling-plus-parent
cluster that `strict` keeps. No-op on this genome; carried into the release note
and pinned by a unit test (§3.5).

**Q3 — confirm the 42 are all chains.** Implied by 459 = 355 + 42 + 62, and
worth an explicit check when re-deriving, because it is exactly the population
Q2 is about.

## 8. Acceptance checklist

- [ ] `--annotation_conflict strict` byte-identical to 0.6.0.0 on all test data.
- [ ] Sorted-order tie-break preserved, covered by a deliberate two-way tie.
- [ ] `nested` on the request's genome: 459 kept (355 + 42 + 62), 5 dropped,
      +62 sequences / +68,162 bp.
- [ ] 69 clusters promoted; 50 to `chromovirus`, 19 to `Chlamyvir`; 0 / 10 / 59
      by prior outcome.
- [ ] Every promoted label is a descendant of that cluster's LCA (asserted).
- [ ] Promotion never fires with two sibling lineages present (asserted).
- [ ] Sibling-plus-parent cluster dropped under `nested` (§3.5), pinned by test.
- [ ] README and changelog state that `nested` can remove sequences, not only add.
- [ ] Cross-superfamily mixes still dropped under `nested`.
- [ ] R3's element collapse covered by unit test — the genome does not exercise it.
- [ ] Summary line printed.
- [ ] `tests/library.sh` in CI.
