# MSA rescue (tier 3) — implementation plan

**Spec / data backing:** the MSA-rescue validation runs in §4.5 of
`docs/refine_v2_analysis.md`-adjacent notes (`tmp/msa_validation_at.tsv`,
`tmp/msa_validation_g2.tsv`).  Headline numbers replicated below for
reference.
**Scope:** add a third refinement tier after parasail-inner and the
existing MAFFT cluster-level fallback — a **per-side MSA rescue** that
fires only on sides the prior tiers left unrefined.  TSD-loss revert
rule already in place gates the result.
**Testing discipline:** each step is a separate commit; smoke +
`./tests.sh refine` must pass locally before each commit.

---

## 1. Files touched

```
utils/
├── refine_boundaries.py         [moderate edits — new tier orchestration]
├── refine_mafft_fallback.R      [tiny edit — widen --snap_window default
│                                 from 5 to 20 bp]
└── (no other utils touched)

dante_ltr_refine                  [unchanged]

tests/
└── refine.sh                    [add tier-3 smoke assertions]

docs/
├── refine_msa_rescue_plan.md    [this file]
└── refine_v2_analysis.md        [append a §9 with the MSA validation
                                  numbers; no rewrite needed]

version.py                       [no bump in this work — comprehensive
                                  evaluation must complete first per
                                  user direction]
```

No changes to `parasail_boundary.py`, `tsd_redetect.py`, or any
downstream tool (`build_ltr_library.R`, `dante_ltr_solo`).

---

## 2. Algorithm changes (semantic summary)

The refinement pipeline gains a third tier that operates **after** the
existing MAFFT cluster-level fallback:

```
tier 1   parasail inner-pool (current default; sets final_corrected_g
         when inner_motif_ok)
tier 2   MAFFT cluster-level fallback (existing — fires on clusters
         with low inner-pool validation rate; sets final_corrected_g
         on members the fallback can resolve)
tier 3   *** NEW *** per-side MSA rescue
         - runs on EVERY qualifying cluster (not just low-validation)
         - applies only to sides where final_corrected_g is still None
         - reuses the same per-cluster MSA the tier-2 fallback already
           builds; we just don't gate on the cluster validation rate
         - the snap-to-TG/CA window widens from ±5 to ±20 bp because
           the MSA change-point typically lands 5–15 bp inside the
           true terminus (validated empirically — see §3 numbers)
tsd-gate per-element TSD recheck + revert rule (unchanged; already
         catches any tier-3 proposal that destroys an existing TSD)
```

Key semantic point: the tiers are **additive**, not exclusive.  Tier 1
is the dominant contributor (handles ~90 % of sides on g2); tier 2
catches a few clusters that fail tier 1 broadly; tier 3 finishes
individual sides that the prior tiers left empty.  The TSD gate is the
single safety net.

**New confidence label `msa_rescue`** for sides whose `final_corrected_g`
came from tier 3.  Label semantics:

| label | source |
|---|---|
| `dual` | tier 1, inner+outer pools agree ≤ 5 bp on coord |
| `divergent` | tier 1, both pools have motif but disagree > 5 bp |
| `inner_only` | tier 1, outer pool failed but inner did |
| `mafft` | tier 2, cluster MAFFT fallback resolved this side |
| **`msa_rescue`** | tier 3, per-side MSA rescue resolved this side |
| `unrefined` | no tier produced a valid coord |

(The `mafft` label is added/separated from the existing
`Refinement_Method=mafft` so the confidence column is no longer used to
distinguish source — see §6 for the corresponding GFF3 schema update.)

**Optional confidence-enhancer flag** (deferred from §10 of this plan):
an `MSA_Agree` per-side GFF3 attribute that compares the MSA call (when
available, regardless of which tier was used) against `final_corrected_g`
within ±5 bp.  Read-only label, no behaviour change.  Decision: include
in this work.

---

## 3. Validation numbers — what tier 3 buys

(From `tmp/msa_validation_*.tsv`; full tables in
`docs/refine_v2_analysis.md` §9 to be added.)

| | at (1 413 TEs in qual. clusters) | g2 (2 793 TEs) |
|---|---:|---:|
| MSA call within ±20 bp of TG/CA | 99.3 % | 99.1 % |
| MSA agrees with v2 inner-pool (`dual`) ≤ 5 bp | 85.5 % | 74.1 % |
| MSA agrees with v2 inner-pool (`divergent`) ≤ 5 bp | 88.1 % | 71.2 % |
| Rescue candidates (sides MSA can resolve where tier 1+2 failed) | 117 | 362 |
| **Net rescues after TSD gate** | **85** | **276** |
| of which gain a TSD that DANTE_LTR didn't find | 30 | 61 |
| TSD-loss attempts blocked by gate | 25 | 86 |

Cost: ~4–10 min wall-clock for full-dataset MSA on the at/g2 scale.

---

## 4. Implementation order

Three commits, each leaving the tool runnable + tests passing.

### Step 1 — Widen the R snap window default

`utils/refine_mafft_fallback.R` currently snaps within ±5 bp.  The MSA
change-point typically lands a few bp inside the actual TG/CA, so the
±5 fails the snap on ~10 % of calls that ARE valid within ±20.

1.1 Change `snap_window <- 5L` (line ~250) to `snap_window <- 20L`.
    Add comment: *MSA change-point lands inside the conserved region
    by 5–15 bp typically; ±20 covers both ends of that range without
    producing false motif matches outside the LTR.*
1.2 No CLI change to the R script (snap window not currently exposed
    as a flag — keep that simple).
1.3 Sanity: re-run `tests/refine.sh`.  Tier-2 MAFFT fallback should
    now produce a slightly higher motif-OK rate when it does fire,
    but the fixture is small enough that the smoke assertions hold
    either way.

**Commit 1:** `refine: widen MAFFT snap window 5 → 20 bp to match
empirical change-point offset`.

### Step 2 — Tier 3 wiring in `refine_boundaries.py`

The substantive change.  All edits are inside one Python module.

2.1 **`RefRecord` extension.**

```python
# tier 3 rescue diagnostics
msa_corrected_g: Optional[int] = None
msa_motif_ok: Optional[bool] = None
msa_agree_with_final: Optional[bool] = None  # |msa - final| <= 5 (or None
                                              # when either is missing)
```

The existing `mafft_corrected_g` / `mafft_motif_ok` fields stay, and
they specifically refer to **tier 2** outputs (cluster-level fallback).
The new `msa_corrected_g` / `msa_motif_ok` are **tier 3** (per-side
rescue and per-side MSA reading).

2.2 **`run_msa_per_cluster(records, cluster_members, …)`.**  New
function, structurally similar to `run_mafft_fallback()` but:

- Always runs on a cluster regardless of validation rate.
- Reuses the same R script (`utils/refine_mafft_fallback.R`) — the R
  script already produces one row per member with `corrected_g` +
  `motif_ok`.
- Rescue rule: for each record where `final_corrected_g is None` AND
  the R script returned a valid motif-OK coord, set
  `final_corrected_g`, `refinement_method = "msa_rescue"`,
  `confidence = "msa_rescue"`.
- Always populates `msa_corrected_g` / `msa_motif_ok` on every
  record (independent of whether rescue fired) so the
  enhance-confidence flag (`MSA_Agree`) can be computed downstream.

2.3 **Driver wiring in `refine_all()`.**

```python
# After tier 1 (parasail per-cluster) and tier 2 (run_mafft_fallback)
if args.msa_rescue:
    for cid, recs in cluster_records.items():
        run_msa_per_cluster(recs, cl_mems_by_id[cid], cid, …)

# Then apply the existing TSD gate on the combined record set.
apply_tsd_gate_per_element(all_records, te_meta, genome, genome_lens, …)
```

Tier 2 may already have populated MSA outputs for the clusters it ran
on; tier 3 should *not* duplicate the R subprocess for those.  Keep a
per-cluster cache: if `cluster_index[cid]["mafft_invoked"]` is True,
parse the same R output for the rescue step.  (Practically: refactor
`run_mafft_fallback` and `run_msa_per_cluster` to share an inner
helper that runs MAFFT + parses output, returning both the records
update and the per-member dict.)

2.4 **`MSA_Agree` computation.**  After all tiers + TSD gate run, for
each record with both `msa_corrected_g` and `final_corrected_g`
populated, set `msa_agree_with_final = abs(msa - final) <= 5`.

2.5 **Confidence aggregation on TE row** (`emit_refined_gff3`).
Update `_CONF_RANK`:

```python
_CONF_RANK = {
    "unrefined":  0,
    "msa_rescue": 1,
    "inner_only": 2,
    "divergent":  3,
    "dual":       4,
    "mafft":      2,   # same tier as msa_rescue / inner_only
}
```

(Adding `mafft` and `msa_rescue` at the same rank as `inner_only`,
since all three are "single-signal-validates" sources.)

2.6 **TE-row method aggregation** stays the same — any side with
non-`none` method "wins" the TE-row method label.

**Commit 2:** `refine: tier-3 MSA rescue for sides unrefined by
parasail+MAFFT cluster fallback`.

### Step 3 — GFF3 / TSV / JSON schema updates + CLI

3.1 **GFF3 schema** (per-LTR row and TE row):

Existing per-LTR attributes gain:
```
MSA_g=<int>                         # always populated when MSA ran
MSA_Motif_OK=TRUE|FALSE|NA          # always populated when MSA ran
MSA_Agree=TRUE|FALSE|NA             # |msa - final| <= 5
```

The TE row carries 5'/3' comma-separated forms:
```
MSA_g=<int_5>,<int_3>
MSA_Motif_OK=<b_5>,<b_3>
MSA_Agree=<b_5>,<b_3>
```

`Refinement_Method` and `Refinement_Confidence` keyword sets gain
`msa_rescue` (and the existing `mafft` is unchanged).

3.2 **Per-element TSV** gains columns (after the existing tier-1 +
tier-2 ones):

```
msa_corrected_g, msa_motif_ok, msa_agree_with_final
```

3.3 **Cluster manifest** gains:

```
msa_invoked          (TRUE iff MSA was run on this cluster)
msa_rescue_count     (count of sides where confidence == "msa_rescue")
```

3.4 **`run.json`** gains under `counts`:

```
n_refined_msa_rescue
n_msa_calls_total              (count of sides with msa_corrected_g != None)
n_msa_agree_with_final         (count of sides where MSA_Agree=TRUE)
```

3.5 **CLI** (`utils/refine_boundaries.py:build_arg_parser`):

```
ap.add_argument("--no_msa_rescue", action="store_true",
                help="Disable tier-3 per-side MSA rescue (default: enabled).")
ap.add_argument("--msa_snap_window", type=int, default=20,
                help="Snap window applied to MSA-derived coords (bp; default 20).")
```

The `--msa_snap_window` is plumbed to the R script (currently the
script has it hard-coded; expose as a CLI arg in the same Step 1 that
widens the default — this lets us tune from Python without re-editing
R).

**Commit 3:** `refine: GFF/TSV schema + CLI for MSA rescue, with
MSA_Agree confidence enhancer`.

---

## 5. Testing

5.1 **Smoke (`tests/refine.sh`).**  Add assertions to the parasail-only
stage (1a):

- `Refinement_Method=msa_rescue` may appear (not required since the
  fixture is small and tier 1 may resolve everything; just check the
  schema for `MSA_g=` attribute presence).
- `MSA_Agree=` attribute appears on at least one row.
- `--no_msa_rescue` produces a refined GFF3 with **no** `msa_rescue`
  in `Refinement_Method` (negative assertion).

5.2 **Acceptance criteria on at + g2** (run after Commit 3):

| dataset | predicted gain |
|---|---|
| at | ~85 net rescues (from 0 → 85 sides with `msa_rescue` confidence) |
| g2 | ~276 net rescues |
| both | TSD-loss outcome must remain 0 (gate intact) |

5.3 **Performance.**  Tier 3 adds a Rscript subprocess per qualifying
cluster.  On at this is +~4 min over tier 1 (8 min total wall-clock);
on g2 about +9 min.  Acceptable.  Document in commit message.

---

## 6. Migration / compatibility notes

- The new `msa_rescue` keyword needs to be accepted by
  `utils/build_ltr_library.R`'s "validated" filter (currently
  `dual / divergent / inner_only / high / medium`).  Add `msa_rescue`
  and `mafft` to the accepted set.  One-line edit; bundled with
  Commit 2.
- The new GFF3 attributes (`MSA_g`, `MSA_Motif_OK`, `MSA_Agree`) are
  appended; downstream readers that don't recognise them ignore them.
- `--no_msa_rescue` keeps the existing v2 behaviour byte-for-byte
  (subject to the snap-window widening from Commit 1, which only
  affects tier 2 coords for clusters where it fires).

---

## 7. Risks

1. **Tier 3 + tier 2 redundancy.**  When tier 2 fires on a cluster, the
   per-side MSA outputs are already produced.  Step 2.3 must NOT call
   the R subprocess twice.  Mitigation: shared inner helper.
2. **Performance on very large clusters.**  MSA is O(N · L) per
   cluster.  On clusters > 500 LTRs, runtime per cluster could exceed
   the existing `--mafft_timeout` (default 600 s).  Mitigation: keep
   the per-cluster timeout; MSA-failure on a giant cluster just means
   no rescue for those members.
3. **`msa_agree_with_final` semantics for divergent / dual.**  When
   final = inner_pool's coord and MSA disagrees, `MSA_Agree=FALSE`
   even though the inner-primary policy chose correctly (data shows
   inner is right ~85 % of the time on disagreements).  Document the
   attribute as "MSA's independent vote, not ground truth" so
   downstream consumers don't misuse it as a "this is wrong" signal.
4. **Snap window widening (5 → 20).**  Affects tier 2 outputs as well.
   Validation showed tier 2 was already barely firing; the change only
   helps the few clusters where it does fire.  Low risk.

---

## 8. Open follow-ups (after this work lands)

- **Confidence enhancer for `divergent`.**  The MSA_Agree attribute
  gives data; once we see how often divergent + MSA-agrees-with-inner
  produces a clean refinement, consider promoting that subset to a
  higher confidence label.  Deferred — needs comprehensive evaluation
  numbers first.
- **Joint-MSA refinement (option C-3 from the v2 plan).**  Still
  deferred to v3.  Tier 3 here is a per-side rescue; joint-MSA would
  be element-level coupled refinement and is a bigger algorithmic
  change.

---

## 9. Acceptance gate (before any version bump)

Per user direction: no version bump in this work.  The comprehensive
evaluation that triggers the next version bump should include tier 3
results across at + g2 + at-fixture, and should at minimum confirm:

- TSD-loss == 0 on both datasets after tier 3 (i.e. gate still
  intact).
- Net rescue counts match the pre-implementation predictions
  (~85 on at, ~276 on g2).
- `dante_ltr_to_library --refined_gff3` and `dante_ltr_solo
  --refined_gff3` round-trip the new schema.
