# Core-domain mode (`--mode core`) — design document

**Status:** design draft, no code changes yet.
**Related:** `docs/fallback_classification_design.md` (a different, lighter
answer to the same problem — see §11).

---

## 1. Problem

`dante_ltr`'s current detector is *lineage-keyed*. A domain cluster is
only considered a candidate element when

1. all its domains carry the **same** `Final_Classification`, and
2. that classification is a **row of `lineage_domain_order.csv`** (one of
   37 REXdb lineages), and
3. the observed domain order matches that lineage's full 5–6 domain
   order within `max_missing_domains`.

On genomes far from REXdb's reference set, DANTE frequently cannot
resolve a domain past `Class_I/LTR/Ty3/gypsy` or
`…/gypsy/chromovirus`, and neighbouring domains of one real element
often get *different* best hits. Both facts break requirement (1) and
(2), so the lineage-keyed clustering never forms, and complete elements
that are structurally obvious are reported only as rank `D`.

The `--fallback_mode` flag (already implemented) softens this by
demoting every classification to a coarse depth and supplying a matching
coarse constraints table. That helps, but it keeps the same structural
requirement: **the full GAG–PROT–…–INT domain complement must be
present and ordered**. Elements whose GAG or PROT is too diverged for
DANTE to call are still missed, and the coarse table still forces one
canonical accessory-domain order per superfamily.

## 2. The idea

Drop the accessory domains from the *detection* criterion entirely.

Use only the three **core** domains — **RT, RH, INT** — whose
*relative order* is the defining structural difference between the two
LTR superfamilies:

```
element 5' ────────────────────────────────────────────── 3'
Ty3/gypsy :  [LTR]  GAG PROT   RT  RH  (aRH)  INT  (CHD)  [LTR]
Ty1/copia :  [LTR]  GAG PROT   INT  RT  RH                [LTR]
                              └──── core ────┘
```

The ordered core triplet is the **seed**. Everything else — GAG, PROT,
aRH, CHD, CHDCR — is demoted from *evidence for detection* to

- **search-space material** — the LTR must lie beyond them, so they
  must be traversed, not stopped at, and
- **classification material** — after the element is delimited, they are
  re-read to refine the call (§7).

Consequences:

- Detection no longer needs lineage-depth classification at all. The
  superfamily falls out of the domain *order*, which is structural and
  therefore robust across large phylogenetic distance.
- A gypsy element whose GAG is undetectable is still found, as long as
  RT, RH and INT are.
- The search space for the LTR is larger than in lineage mode, because
  it must span the unknown accessory region as well as the UTR — but
  measurement (§5.2) shows only about half as large as this document's
  first estimate assumed. It is affordable because we search for the
  **closest** direct repeat pair and gate on TG…CA + TSD (§6.4).
- The price is specificity: we must replace "the domain complement
  matched a known lineage" with explicit structural gates. §6.4 and §9
  list them.
- The other price is that the core must be **complete**. Requiring all
  three domains is what makes the superfamily call sound, but 2–25 % of
  validated elements do not carry all three through the domain filter,
  depending on the genome (§5.2). That ceiling is set by
  `dante_filtering`, not by the clustering this mode replaces, and
  §6.0.2 relaxes it for core domains specifically.

The central premise — that the core order is invariant and determines
the superfamily — is not assumed. It was measured on 51 875 validated
elements across three genomes and held in every one (§5.2).

And the payoff is measured too: on *Draparnaldia*, which REXdb does not
cover, lineage mode returns **zero** complete elements while core-mode
seeding finds **1 984** ordered cores on the same input (§5.4).

**Non-goal:** replacing lineage mode. On REXdb-covered species, lineage
mode's stricter criterion is the better one and stays the default.

---

## 3. User-visible surface

### 3.1 New flag

```
dante_ltr -g DANTE.gff3 -s genome.fasta -o out --mode core

--mode {lineage,core}   default: lineage
    lineage : current behaviour — lineage-keyed full domain complement
    core    : seed on the ordered RT/RH/INT core, classify afterwards
```

Core-mode tuning flags. `[table]` means the default is the
per-superfamily value from `databases/core_domain_order.csv` (§5); a CLI
value overrides the table for both superfamilies:

```
--core_max_gap N        max bp between consecutive core domains   [table]
--core_max_span N       max bp from first to last core domain     [table]
--min_ltr_length N      minimum LTR length to accept              [100]
--max_te_length N       maximum element length                    [35000]
--core_require_tsd      reject elements without a TSD  (rank DL)  [off]
--min_relative_length_core F   relaxed Relat_Length for RT/RH/INT  [0.3]
--min_similarity F      dante_filtering similarity threshold      [0.4]
--split_max_gap N       max bp between fragments of one domain    [50]
```

`--min_relative_length_core` is the core-mode sensitivity knob: the
ordered triplet is strong enough joint evidence to admit individually
weak core domains (§6.0.2). `--min_similarity` exposes what is today a
hard-coded constant in `dante_filtering()` (`ltr_utils.R:951`); its
default is unchanged and §5.3 shows it is non-binding in practice, so it
exists for completeness rather than tuning (O1).

### 3.2 Everything else is unchanged

Same inputs, same output file set, same feature types, same ranks
(`DLTP`/`DLT`/`DLP`/`DL`/`D`), same statistics CSV shape. Downstream
tools (`dante_ltr_summary`, `dante_ltr_to_library`, `clean_ltr.R`,
`dante_ltr_solo`) consume `Rank`, `Final_Classification` and the
feature types; none of them key on `source` or on the new attributes, so
they work on core-mode output without modification.

Core-mode elements are still distinguishable: `source=dante_ltr_core`
in column 2, plus the attributes of §8.

### 3.3 Relation to `--fallback_mode`

Independent and combinable. `--fallback_mode` demotes the *input*
classifications; core mode does not depend on input classification depth
for detection, so combining them is legal but usually pointless. When
both are given, core mode's own classification logic (§7) governs the
output, and a note is printed. Document the overlap in the README;
do not deprecate `--fallback_mode` (it remains the cheap option when
the full domain complement *is* present but under-resolved).

---

## 4. Where the code goes

```
dante_ltr                      (python)  + --mode dispatch, arg plumbing
└─ utils/detect_putative_ltr.R           unchanged
└─ utils/detect_core_ltr.R     (new)     core-mode driver, mirrors the
                                         structure of detect_putative_ltr.R
└─ utils/core_ltr_utils.R      (new)     seeding, search-space, reannotation
└─ utils/ltr_utils.R                     unchanged; sourced by both
databases/core_domain_order.csv (new)    superfamily-level constraints
```

Python side: `--mode` only selects the R script name in
`_detect_ltr_cmd()` (line ~274) and in the single-chunk path
(line ~1591), and appends the core-mode options. **Chunking, the
memory-gated process pool, coordinate remapping, statistics merging,
`add_version_to_gff3` and the fallback preprocessing all stay exactly as
they are** — core mode is a different per-chunk detector behind the same
harness, so it inherits chunk parallelism and large-genome support for
free.

Reused verbatim from `ltr_utils.R` (no edits):

| function | used for |
|---|---|
| `CHD_CHDCR_correction` / `revert_CHDCR_correction` | CHD naming |
| `gff_cleanup_overlaps` | overlapping DANTE hits |
| `dante_filtering` | domain quality filter |
| `add_coordinates_of_closest_neighbor` | neighbour bookkeeping |
| `blast`, `trim2TGAC`, `mask_tandem_repeats` | LTR pair search |
| `evaluate_ltr`, `get_best_ltr`, `get_TE` | LTR + TSD calling |
| `get_te_gff3` | element → GFF3 features |
| `add_pbs`, `add_pbs_hemi` | PBS / tRNA |
| `get_te_rank` | DLTP/DLT/DLP/DL |
| `get_te_statistics`, `get_te_sequences` | outputs |
| `trim_gr`, `merge_gr`, `add_info_about_flanking_sequences` | assembly of final GFF3 |

New code is confined to §6.0.1–§6.3 and §7. That is the whole point of
this design: the expensive, well-tested parts (LTR/TSD/PBS) are
untouched, and only the *candidate generation* and *classification*
stages are replaced.

Two measurement tools support the design and are not part of the
shipped pipeline:

| script | input | answers |
|---|---|---|
| `utils/calibrate_core_constraints.py` | `dante_ltr` output GFF3 | the constraints table (§5.1) |
| `utils/measure_core_seeding.py` | raw DANTE GFF3 | fragmentation and filter thresholds (§5.3) |

---

## 5. `databases/core_domain_order.csv`

Same tab-separated shape as `lineage_domain_order.csv`, keyed by
superfamily instead of lineage, with the core-specific columns added.
User-overridable via `--te_constrains` (which in core mode is read as
this table, not the lineage one).

Values below are **measured**, not guessed — see §5.1 for the procedure
and §5.2 for the results. They are still preliminary: three genomes,
pending the Darwin runs.

```
Superfamily	core_order	offset5prime	offset3prime	core_max_gap	core_max_span	ltr_length	accessory_5	accessory_3
Class_I/LTR/Ty1_copia	INT RT RH	9000	8500	2000	6000	100	PROT GAG	-
Class_I/LTR/Ty3_gypsy	RT RH INT	17000	14000	2000	6000	100	PROT GAG	aRH CHD CHDCR
```

Notes on the columns:

- `offset5prime` / `offset3prime` are distances **from the core**, not
  from the first/last domain of the full element as in
  `lineage_domain_order.csv`. They are the largest value observed across
  the three measured genomes, rounded up — deliberately not
  quantile-trimmed, for the reason given in §5.2.
- `core_max_gap` and `core_max_span` are set at roughly twice the q99
  instead, because an over-large gap cap can join domains from adjacent
  elements. The gap cap comfortably covers RH→INT across an aRH in
  Tat-like elements, which is the largest legitimate core gap.
- `ltr_length` is the *expected* LTR length fed to `blast()`, which
  filters alignments shorter than `0.8 ×` it. It is deliberately **not**
  taken from the measurement (see §5.1) and is pinned at the 100 bp
  floor, matching the coarse fallback tables.
- `accessory_5` / `accessory_3` list which non-core LTR-RT domains may
  legitimately appear on each side of the core **in element
  orientation**; this drives the transparency rule of §6.3. The copia
  3' side is empty — RH is the last domain of a copia element — which
  means a copia core's 3' walk blocks at the first domain it meets.

### 5.1 Calibrating the table from annotated genomes

The quantities the table needs are all directly observable in existing
`dante_ltr` output, because a rank `DLT`/`DLTP` element has both its
LTR boundaries and its domain positions in one GFF3. Measuring on
elements that lineage mode already validated gives ground truth for
where the LTR sits relative to the core — exactly the number core mode
has to guess.

`utils/calibrate_core_constraints.py` (written alongside this document;
stdlib-only, streams its input) takes one or more `dante_ltr` output
GFF3s and emits the table. Over the Darwin runs:

```bash
utils/calibrate_core_constraints.py \
    -g '/nfsroot/projects/darwin/runs/*/dante_ltr*.gff3' \
    -o databases/core_domain_order.csv \
    --report core_constraints_report.tsv \
    --margin 1.5
```

Per superfamily, over elements of rank `DLT`/`DLTP` only (TSD present ⇒
boundaries trustworthy), it measures:

| measured quantity | sets |
|---|---|
| gap between consecutive core domains, in element orientation | `core_max_gap` |
| first core domain start → last core domain end | `core_max_span` |
| 5'LTR start → core start | `offset5prime` |
| core end → 3'LTR end | `offset3prime` |
| element length | `--max_te_length` default |
| which accessory domains occur 5' vs 3' of the core, and how often | `accessory_5` / `accessory_3` |

Two measurement caveats, both of which change how the output is used:

- **The offsets are censored.** The only elements available to measure
  are those lineage mode already found, i.e. those whose LTR fell
  inside lineage mode's *own* search window. Elements with unusually
  long UTRs are systematically absent, so the observed distribution is
  a lower bound. Hence the 99th percentile is taken with a **1.5×
  margin** (`--margin`). Erring high is cheap: a wider window costs
  BLAST time, while a window that is too narrow loses the element
  outright.
- **LTR width must not be taken from the measurement.** It would be
  circular — lineage mode found these LTRs *using* the per-lineage
  `ltr_length` prior in `lineage_domain_order.csv`, so the observed
  widths inherit that prior and cannot justify it. The script reports
  the distribution for information, but `ltr_length` stays pinned at
  the 100 bp floor.

The report also answers two questions the design rests on:

- **Is the core order actually invariant?** Count elements whose core
  domain order contradicts their lineage-assigned superfamily. This is
  the central premise of core mode; if it is not ≈0, the design needs
  revisiting before any code is written.
- **How often is the core complete?** Fraction of validated elements
  carrying all three of RT/RH/INT above the filtering threshold. This
  is core mode's sensitivity ceiling on well-covered genomes, and it
  bounds what §10's concordance test may demand.

### 5.2 Measurement on three annotated genomes

Run over the three distinct annotated genomes in `test_data/`
(`g1_dante_ltr.gff3` and `at/Alyr_dante_ltr.gff3` are the same file):

| | g1 / Alyr | g2 | g3 / Pisum |
|---|---:|---:|---:|
| elements at rank DLT/DLTP | 1 483 | 4 642 | 47 529 |
| copia: clean core | 449 (79 %) | 3 685 (98 %) | 9 294 (93 %) |
| gypsy: clean core | 686 (75 %) | 862 (98 %) | 36 899 (98 %) |
| **core order vs superfamily** | **1135 / 1135** | **4547 / 4547** | **46193 / 46193** |

Per-superfamily geometry, q99 (and observed max):

| | | g1 | g2 | g3 |
|---|---|---:|---:|---:|
| copia | core gap | 949 | 1 000 | 828 |
| | core span | 3 283 | 3 103 | 2 931 |
| | offset5′ | 3 335 (4 072) | 3 516 (4 778) | 4 894 (8 832) |
| | offset3′ | 3 040 (3 475) | 2 207 (3 709) | 4 937 (8 490) |
| | TE length | 8 508 | 8 468 | 12 164 |
| gypsy | core gap | 1 008 | 774 | 804 |
| | core span | 3 134 | 3 634 | 2 919 |
| | offset5′ | 5 212 (7 949) | 9 542 (10 016) | 14 164 (16 992) |
| | offset3′ | 7 666 (11 975) | 7 044 (7 559) | 11 590 (13 996) |
| | TE length | 14 779 | 17 945 | 28 190 (29 951) |

**1. The core-order premise holds exactly — 51 875 / 51 875.** Across
three genomes, every element's core order agrees with its lineage
assignment, and not one accessory domain's superfamily disagrees with
its element's. Calling the superfamily from domain order alone is
sound. This is the design's central assumption and it is now measured,
not assumed.

**2. Offsets are strongly genome-dependent — measuring one genome would
have been actively misleading.** Gypsy `offset5prime` q99 runs 5 212
(Alyr) → 9 542 → 14 164 (Pisum), a 2.7× spread tracking element size
(Pisum carries the large Ogre/Tat elements). An earlier draft of this
document measured only Alyr and concluded the offsets could be halved;
across three genomes that is wrong, and the original extrapolation from
`lineage_domain_order.csv` was close to right. The table takes the
**maximum observed value across genomes**, rounded up.

Erring generous is the correct bias here, and cheaply so: `get_TE()`
sorts BLAST hits innermost-first, so widening the window only adds
candidate pairs *outside* the true LTRs, which sort behind the real one.
A window that is too wide costs BLAST time; one that is too narrow loses
the element outright. This is why the offsets are not quantile-trimmed.

**3. Core gap and span are stable across genomes** — gap q99 ≈ 800–1 000
and span q99 ≈ 2 900–3 600 everywhere, unlike the offsets. The rare tail
(gap max 5 583) is set aside deliberately: a large gap cap risks joining
domains from adjacent elements, which is a correctness problem, whereas
a missed outlier is only a sensitivity loss.

**4. Core completeness varies 75 % → 98 %,** and the Alyr figure is the
outlier, not the rule. It still sets a real ceiling on some genomes, so
§10's concordance target stays conditional on the core being present.

**5. `get_best_ltr()`'s hard 30 kb element cap is being hit.** Pisum
gypsy elements reach 29 951 bp with q99 at 28 190 — the distribution is
clipped against the constant in `ltr_utils.R:711`. `--max_te_length`
should default to 35 000 in core mode rather than inheriting 30 000.

Resulting table (maximum across genomes, rounded up; gap and span at
2× q99):

```
Superfamily	core_order	offset5prime	offset3prime	core_max_gap	core_max_span	ltr_length	accessory_5	accessory_3
Class_I/LTR/Ty1_copia	INT RT RH	9000	8500	2000	6000	100	PROT GAG	-
Class_I/LTR/Ty3_gypsy	RT RH INT	17000	14000	2000	6000	100	PROT GAG	aRH CHD CHDCR
```

Caveat: three genomes, all reasonably REXdb-covered. The Darwin runs
should widen this before the table is frozen — the command in §5.1 pools
them, and per-genome runs give the spread that matters for the offsets.

---

### 5.3 Fragmentation and filter thresholds — raw DANTE measurement

§5.2 cannot answer two questions, because both concern annotations that
are absent from `dante_ltr` output by construction (§6.0.1, §6.0.2).
`utils/measure_core_seeding.py` measures them on raw DANTE instead,
applying a configurable filter, optionally merging split annotations,
and enumerating seeds exactly as §6.2 specifies.

Measured on the Pisum JI1006 raw DANTE (1 439 748 domains, 58
sequences); its lineage-mode annotation is `test_data/g3`:

```
domains in file                 : 1 439 748
passing filter (Sim .4/RelL .6) : 1 045 011  (72.6 %)
core domains passing            : RT 173 070   RH 186 402   INT 173 235
ordered core seeds              : 128 010  (gypsy 97 113 / copia 30 897)
seed span q50/q95/q99/max       : 2 858 / 2 916 / 2 919 / 4 776

lineage mode on the same genome : DLTP 37 233  DLT 10 296  DLP 13 009
                                  DL 23 870    D 161 320
```

#### Fragmentation is real, hits RT hardest, and is rarer than expected

Merging must run **before** filtering — each fragment's `Relat_Length`
is reduced by construction, so filtering first destroys exactly the
annotations we are trying to rejoin. With `--split_max_gap 4000` so the
distribution is not clipped by the cap:

```
same-name adjacent pairs : 6 563 classified split
split by domain          : RT 2441, INT 1658, RH 765, TPase 580, GAG 460, ...
weaker fragment Relat_Length (q05/q50/q95) : 0.09 / 0.26 / 0.47
split gap histogram (bp) :
     0-20    244        20-50    129        50-100   136
   100-200   419       200-500  1112       500-1k    818
     1k-2k  2421         2k-4k  1255
```

RT is the most frequently affected domain, as predicted. And the
mechanism is confirmed: **the weaker fragment's `Relat_Length` is 0.26
at the median and 0.47 at q95 — below the 0.6 filter in essentially
every case.** The short half is silently dropped, the long half
survives, and the reference-domain coverage that would have justified
both is never summed.

But the raw count of 6 563 is mostly background, and the histogram shows
it. Converting each bin to a per-bp density exposes two populations:

```
bin        0-20   20-50  50-100  100-200  200-500  500-1k  1k-2k  2k-4k
pairs/bp   12.2    4.3     2.7      4.2      3.7     1.6     2.4    0.63
```

A frameshift or stop codon breaks a domain with a gap of *tens* of bp,
and that is where the excess sits: the 0–20 bin runs at 12.2 pairs/bp
against a ~2–4 pairs/bp background that is roughly flat out to 4 kb.
Background of that shape is what unrelated pairs produce — a wider
window simply catches more of them. Subtracting it leaves **roughly
200–250 genuine frameshift splits genome-wide, against 173 854 RT
domains: about 0.1 %.**

The kb-scale pairs are not frameshifts. They are either domains
interrupted by a nested insertion or, more likely at these distances,
two degraded neighbouring elements whose partial domains happen to tile
the reference. Complementary tiling alone is too weak a test once the
gap exceeds a few hundred bp. **`--split_max_gap` should therefore
default to ~50 bp**, where the excess is concentrated; the apparent
"gains" at larger caps are spurious merges. The measurement at
`--split_max_gap 1000` (+204 seeds) and at 4000 (+892 seeds) should both
be read as mostly background, not sensitivity.

#### So merging is correctness work, not a sensitivity lever

```
seeds                             128 010  ->  128 214   (+0.2 %, cap 1000)
seeds blocked by a core domain        379  ->      389
```

Merging moves seeding by a fraction of a percent, and the failure it is
meant to prevent — a leftover RT half blocking the outward walk —
affects 389 of 128 214 seeds (0.3 %). The reason is the same in both
cases: the weak half usually fails the filter, so it neither seeds nor
blocks.

That is a weaker result than the hypothesis suggested, and it should be
stated plainly: on this genome, fragment merging fixes a real,
well-diagnosed failure mode affecting ~0.1–0.3 % of elements. It is
cheap and it is correct, so §6.0.1 keeps it — but it must not be sold as
a sensitivity feature.

The natural next guess — that frameshifted elements would be commoner on
a genome REXdb does not cover — was tested on Drapa (§5.4) and **does
not hold**: there the short-gap excess disappears entirely into
background. Fragment merging therefore has no genome in evidence where
it moves sensitivity.

#### The transparency rule is the load-bearing one

```
seeds with a domain within 500 bp of a seed edge : 62 674 / 128 214 (49 %)
  PROT 44 267   CHD 17 289        <- transparent by design (§6.3)
  GAG 495   CHDCR 145   aRH 65    <- transparent by design
  RH 166   RT 123   INT 100       <- blocking (389 total, 0.3 %)
```

**Half of all seeds have a domain within 500 bp of the core**, and 99 %
of those are PROT/CHD/GAG — the accessory domains §6.3 traverses. Had
the design kept lineage mode's naive "stop at the nearest neighbour"
clamp, half of all core seeds would have had their LTR window truncated
to nothing. This is the strongest quantitative justification for §6.3
in the document.

#### Similarity is the wrong knob; Relat_Length is the lever

```
Similarity    domains passing    seeds
      0.25          1 046 080    128 214
      0.40          1 046 079    128 214        <- no effect at all

Relat_Length  domains passing    seeds         vs 0.6
      0.30          1 129 131    141 573       +10.4 %
      0.40          1 092 473    133 784        +4.3 %
      0.50          1 068 893    130 957        +2.1 %
      0.60          1 046 079    128 214            —
```

Similarity between 0.25 and 0.40 changes the seed count by **zero** on
this genome — it is entirely non-binding. `Relat_Length` is what
actually gates core-domain admission, and relaxing it to 0.3 yields
10 % more seeds.

This corrects the framing of O1 and of §6.0.2, both of which were
written around `min_similarity`. The core-relaxed threshold should be
**`--min_relative_length_core`**; `--min_similarity` should still be
exposed for completeness, but nothing suggests moving it. §5.4 confirms
both findings on a genome at the opposite end of the coverage range.

#### Caveats

- **Seeds are not elements.** A seed becomes an element only if an LTR
  pair is found (§6.4). These counts bound core-mode sensitivity from
  above; they do not predict output.
- **The "duplicate pair" count (66 086) is not meaningful as measured.**
  Run on raw DANTE, it conflates rival overlapping hits at one locus —
  which the real pipeline removes in `gff_cleanup_overlaps()` first —
  with genuinely separate domains. Only the *split* classification,
  which requires complementary reference-domain tiling, is reliable here.
- **One genome, and a well-covered one.** §5.4 repeats the measurement
  on Drapa, which REXdb does not cover, and refutes the expectation that
  both effects would be larger there.

---

### 5.4 Drapa — a genome REXdb does not cover at all

*Draparnaldia* (a green alga) is the motivating case for this work and
for `docs/fallback_classification_design.md`. Measured on
`.../DRA_2025_07_30/.../output/DANTE/DANTE.gff3` (61 446 domains) and
the matching `DANTE_LTR.gff3`.

#### Lineage mode finds nothing at all

```
dante_ltr output on Drapa, by rank:   D 201
                                      DL 0   DLT 0   DLP 0   DLTP 0
```

**Not one complete element.** `calibrate_core_constraints.py` cannot run
on this genome — there is no rank `DLT`/`DLTP` element to measure — and
the whole `DANTE_LTR.gff3` is 586 KB of rank-`D` domain clusters. This
is the clearest possible statement of the problem core mode addresses,
and it is worth quoting in the README.

#### Core mode has ~2 000 candidate seeds on the same input

```
domains in file                 : 61 446
passing filter (Sim .4/RelL .6) : 24 836  (40.4 %)   [Pisum: 72.6 %]
core domains passing            : RT 3 206  RH 3 476  INT 3 481
ordered core seeds              : 1 984   (gypsy 1 779 / copia 205)
seed span q50/q95/q99/max       : 2 457 / 2 974 / 3 369 / 4 137
```

From **0 complete elements to 1 984 ordered core seeds** on identical
input. Seeds are not elements — each still has to yield an LTR pair
(§6.4) — so this is an upper bound, not a prediction. But the premise
that structural seeding survives where lineage-keyed clustering does not
is now demonstrated on the genome it was designed for.

The composition is consistent with an alga: overwhelmingly gypsy
(1 779 vs 205 copia), and chromoviral — CHD and CHDCR together account
for 4 640 of the raw LTR domains against only 489 GAG.

Note the filter bites much harder here: 40 % of domains pass versus
73 % on Pisum, as expected when every hit is to a distant reference.

#### The transparency rule is not merely important here — it is decisive

```
seeds with a domain within 500 bp of a seed edge : 1 967 / 1 984  (99.1 %)
  CHD 1 445   PROT 518        <- transparent by design (§6.3)
  RH 3   RT 1                 <- blocking (4 total, 0.2 %)
```

On Pisum this was 49 %; on Drapa it is **99 %**. Essentially every core
seed has an accessory domain — usually CHD, immediately 3' of the core
in a chromovirus — sitting within 500 bp of the core edge. With lineage
mode's "stop at the nearest neighbour" clamp, virtually every seed on
this genome would have a useless LTR window. §6.3's traversal rule is
what makes core mode work here at all.

#### Two hypotheses from §5.3 are refuted

**O8 predicted both effects would be larger on a distant genome. They
are not.**

*Fragmentation.* Applying the same per-bp density analysis as §5.3:

```
bin        0-20  20-50  50-100  100-200  200-500  500-1k  1k-2k  2k-4k
pairs         7     17      19       48       96     310    156    165
pairs/bp   0.35   0.57    0.38     0.48     0.32    0.62   0.16   0.08
```

There is **no excess at short gaps** — the density is flat across the
whole range, i.e. what remains is background with no detectable
frameshift signal. At the recommended `--split_max_gap 50` the whole
genome yields ~24 split pairs. Merging moves seeding by +36 seeds
(+1.8 %) and changes the count of core-blocked seeds not at all (4 → 4).

So fragment merging stays in the design on the same footing as §5.3 gave
it — cheap, correct, fixes a diagnosable failure — but there is now no
genome in evidence where it is a sensitivity lever, and the design
should not claim otherwise.

*Filter relaxation.* The `Relat_Length` sweep gains **less** on Drapa
than on Pisum:

```
Relat_Length   domains passing   seeds    vs 0.6      (Pisum, vs 0.6)
      0.20             39 205    2 190    +8.4 %
      0.30             32 581    2 153    +6.6 %            +10.4 %
      0.40             28 922    2 106    +4.3 %             +4.3 %
      0.50             26 255    2 054    +1.7 %             +2.1 %
      0.60             25 022    2 020        —                   —

Similarity 0.20 / 0.30 / 0.40  ->  2 020 seeds at every value
```

`Similarity` is again **exactly** non-binding — identical seed counts at
0.20, 0.30 and 0.40. That now holds on two genomes at opposite ends of
the coverage range, which settles O1 conclusively.

The informative number is the ratio: relaxing `Relat_Length` from 0.6 to
0.2 admits **57 % more domains but only 8.4 % more seeds**. The ordered
triplet requirement absorbs nearly all of the extra material — a weak
domain only produces a seed if it happens to lie in the right order, on
the right strand, within `core_max_span` of two others. That is the
"cumulative evidence" argument of §6.0.2, and it is now quantified:
**the relaxation is cheap in specificity precisely because the triplet
constraint, not the per-domain filter, is what carries the weight.**

On a genome where lineage mode returns zero, +170 seeds for no
structural cost is worth having. `--min_relative_length_core` should
default to **0.3**; `--min_similarity` should not be touched.

---

### 5.5 Head-to-head against lineage mode on Pisum

The measurement §13 said was missing: both modes run on the same input,
counting *elements* rather than seeds. Input is
`test_data/sample_genome.fasta` + `sample_DANTE.gff3` — 87 Mb of Pisum
across 50 contigs, well covered by REXdb (Ogre, SIRE, Tekay, Ivana,
Tork; only 116 of 11 061 LTR domains stop at superfamily depth).
Reproduce with `utils/compare_detection_modes.py`.

```
                    D      DL    DLT    DLP   DLTP    total    runtime
  lineage        2473     110     24     55     54     2716        65 s
  core           2194     251     45    111     93     2694       191 s
```

#### Concordance — the §10 gate passes

```
lineage complete elements (rank > D)          243
core    complete elements                     500
matched within +/-20 bp                       233   = 95.9 %
  ...of which exact, 0 bp on both ends        233   = 100 %
superfamily agreement on matched pairs     233/233
classification depth identical             233/233
```

**95.9 % recovery, and every matched element agrees to the base.** Not
one matched pair differs by even a single bp, so the ±20 bp tolerance in
the §10 gate is never exercised. The order-derived superfamily never
contradicts lineage mode, now over 233 more elements. And core mode
demoted nothing here — on a genome where the evidence supports lineage
depth, the LCA rule reports lineage depth.

Every one of the 243 lineage elements has a complete filtered core on
this genome, so the conditional and unconditional figures coincide; the
2–25 % ceiling of §5.2 does not bite here.

#### Sensitivity — core mode roughly doubles the yield

267 complete elements that lineage mode does not report. They are not
noise: their profile tracks the matched set on every axis, just slightly
more diverged, which is what an element lineage mode misses should look
like.

```
                 n     median length   median LTR identity   TSD    PBS
  matched      233            9 728              90.3 %     32 %   45 %
  core-only    267            9 142              88.0 %     23 %   37 %
```

39 of the core-only elements reach rank `DLTP` and 53 more `DLP`, i.e.
92 carry tRNA/PBS evidence independent of the LTR call. Split evenly
between superfamilies (136 copia / 131 gypsy).

#### The 10 misses are boundary choices, not detection failures

Core mode finds something at 9 of the 10 loci; on 6 it calls a
*different, inner* LTR pair, usually sharing one end exactly:

```
  ctg137:520921-529829   lineage DLP    core 523168-529829  (5' end +2247)
  ctg993:1032128-1049165 lineage DLTP   core 1033012-1049165 (5' end  +884)
```

The cause is window geometry. Lineage mode anchors its search window on
the first domain of the cluster — usually GAG — so the window cannot
reach into the element's own GAG/PROT region. Core mode anchors on the
core, 2–3 kb further in, so that region is inside the window, and
`get_TE()`'s innermost-pair preference can then prefer a repeat found
there over the true LTR.

**Anchoring the window on the outermost traversed accessory domain was
tried and is measurably worse**: it repaired 2 of the 10 and broke 6
others, 95.9 % → 94.2 %. The reason is the symmetric failure — when the
traversed domain actually belongs to the *neighbouring* element, the
anchor jumps past the true LTR and the element is lost entirely rather
than merely mis-bounded. The core-anchored geometry is kept.

A better fix, not attempted here, is to make the preference conditional
rather than the window narrower: keep the wide window, but among
candidate LTR pairs prefer the outermost one still consistent with every
accessory domain the walk traversed. That needs a change to `get_TE()`'s
selection, which is shared with lineage mode, so it belongs in its own
piece of work. **Known limitation: ~2.5 % of elements (6 of 243) get a
5' boundary placed inside the true LTR on a well-covered genome.**

#### Cost

Core mode is ~3× slower (191 s vs 65 s) on this input. That is inherent
rather than a defect: 855 seeds against 365 lineage clusters, each with
a wider BLAST window. Serial per-element domain lookups were also
quadratic in (domains × elements) — replaced with an indexed lookup,
verified byte-identical, worth only ~5 s here but necessary before a
real chunk with 10^5 domains.

---

## 6. Algorithm

### 6.0 Preprocessing

`CHD_CHDCR_correction` → `gff_cleanup_overlaps` → prefilter neighbour
counting → **fragment merging (new, below)** → `dante_filtering` →
`add_coordinates_of_closest_neighbor`.

#### 6.0.1 Merging split domain annotations

A frameshift or an in-frame stop can break one protein domain into two
adjacent DANTE annotations — a left part and a right part, each with a
reduced `Relat_Length`, together tiling the reference domain. Core mode
must collapse these into one logical domain **before** seeding, for two
independent reasons:

- **Seeding.** A split RT presents as `RT RT RH INT`. The triplet is
  still enumerable, but the leftover RT is then a core domain outside
  the seed — which §6.3 classifies as *blocking*, truncating the LTR
  search window to nothing. The element is found and then immediately
  lost.
- **The blocking walk.** Transparency rule 4 (§6.3) refuses a second
  copy of an accessory type, on the grounds that it belongs to the
  neighbouring element. A split GAG trips exactly that rule and blocks
  the walk inside the element it belongs to.

Detection uses the signature that distinguishes a split from a genuine
duplication: the two annotations are adjacent on the same strand, and
their `Best_Hit_DB_Pos` intervals **tile complementary parts of the same
reference domain** in element order. Two real copies instead hit
overlapping parts of the reference. Merged domains span both fragments,
sum their `Relat_Length`, and record `Nfragments` in the output.

`--split_max_gap` defaults to **50 bp** and `--split_max_db_overlap` to
0.3. The gap default is tight on purpose: §5.3 measures the frameshift
excess as concentrated below ~20 bp, with complementary tiling losing
specificity beyond a few hundred bp, where unrelated degraded neighbours
start to satisfy it.

This is not only a core-mode concern — it is a class of element lineage
mode **cannot find at all**. `clean_domain_clusters()`
(`ltr_utils.R:206`) requires `N_unique_domains == N_domains`, so any
cluster containing a repeated domain name is discarded outright. Every
element with a frameshifted domain is therefore absent from lineage-mode
output, which is also why the §5.2 calibration reports zero repeated
domains: the phenomenon is censored out of the data it measures.
Quantifying it requires raw DANTE input — see §5.3.

#### 6.0.2 Two filter thresholds, not one

The domain filter is load-bearing in two opposite directions, and a
single threshold cannot serve both:

- For **seeding**, a lower threshold is justified. The ordered triplet
  is itself strong joint evidence: a marginal RT only produces a seed if
  it happens to lie in the correct order, on the correct strand, within
  `core_max_span` of an RH and an INT. Three weak-but-consistent
  domains carry more weight than any one of them alone, so requiring
  each individually to clear the single-domain threshold discards real
  elements — precisely the distant-genome case core mode exists for.
- For **blocking**, the threshold must stay high. §6.3 reads survival
  of the filter as proof that a domain is real and therefore delimits
  the element. Admitting marginal domains as blockers would truncate
  windows on noise.

The §5.3 sweep identifies **which** threshold to relax, and it is not
the one O1 asked about. On Pisum, moving `Similarity` from 0.40 to 0.25
changes the seed count by zero — it is entirely non-binding.
`Relat_Length` is what actually gates core-domain admission: 0.6 → 0.3
yields 10 % more seeds.

So the relaxed threshold is **`--min_relative_length_core`** (default
**0.3**, confirmed on both Pisum and Drapa), applied to RT/RH/INT when
enumerating seed candidates;
`--min_relative_length` (default 0.6) continues to govern everything
else, including the blocking walk and reannotation. `--min_similarity`
stays at 0.4 and is exposed for completeness only (O1).

This also compounds with §6.0.1: a fragmented domain's parts each have a
low `Relat_Length` by construction, so merging and this relaxation
address the same failure from two sides — merging restores the summed
coverage, and the lower threshold admits what is still short after
merging.

### 6.1 Core-domain selection

A filtered DANTE domain is a **core candidate** when

- `Name ∈ {RT, RH, INT}`, **and**
- it has LTR-retrotransposon support: either
  `Final_Classification` starts with `Class_I|LTR`, **or** at least one
  entry of `Region_Hits_Classifications` starts with `RT|Class_I|LTR`
  (resp. `RH|…`, `INT|…`).

The second clause matters for the target use case: on distant species a
genuine LTR RT can win its best hit against `Class_I|pararetrovirus` or
a LINE RT while still having LTR hits in its candidate list. The
candidate list is exactly the information `Region_Hits_Classifications`
already carries (it is currently discarded at the top of the `repeat{}`
block in `detect_putative_ltr.R`). Domains accepted only via the second
clause get `Domain_LTR_Support=secondary` and are counted separately in
the log.

All other filtered domains are kept in a second track used for
search-space delimitation (§6.3) and reannotation (§7) — they are never
seed material.

### 6.2 Seeding: finding ordered core triplets

Per sequence, per strand, over core candidates sorted by coordinate:

1. **Runs.** Split into maximal runs where consecutive core candidates
   are on the same strand and separated by ≤ `core_max_gap`.
2. **Pattern enumeration.** Within a run, enumerate index triples
   `i < j < k` such that the domain names, read in **element
   orientation** (coordinate order on `+`, reversed on `-`), spell
   either

   - `RT, RH, INT` → **Ty3/gypsy**, or
   - `INT, RT, RH` → **Ty1/copia**,

   and such that `end(k) − start(i) ≤ core_max_span` and each
   consecutive gap ≤ `core_max_gap`.
3. **Scoring.** Prefer the most compact triple: sort by span ascending,
   breaking ties on `sum(Similarity)` (O2, resolved). The measured core
   span is tight — q99 ≈ 2.9–3.6 kb across four genomes (§5.2, §5.4) —
   so a stretched triple is more likely to be domains borrowed from two
   neighbouring elements than a real core. Span ties are rare in
   practice, so the tiebreak seldom decides anything; the point of the
   rule is that it has **no free parameter to calibrate**, unlike a
   bitscore-minus-λ·span score, which would need tuning per genome —
   exactly the kind of tuning core mode exists to avoid.
4. **Greedy non-overlapping selection.** Take triples in descending
   score, accept one if it shares no domain with an already accepted
   one. This resolves tandem arrays (`RT RH INT RT RH INT` → two
   seeds) and nested/fragmentary cases without special-casing.
5. A run yielding no valid triple produces no seed; its domains remain
   available to the `TE_partial` (rank `D`) track, exactly as in
   lineage mode.

Each seed carries: `seqnames`, `strand`, core span, the three domain
rows, and `Superfamily` — derived **purely from the order**, with no
reference to any `Final_Classification`.

### 6.3 Search-space delimitation

This is the part that cannot reuse lineage mode. In
`detect_putative_ltr.R`, `get_ranges_left()`/`get_ranges_right()` clamp
the BLAST window at the nearest neighbouring domain, which works there
because the element's own GAG/PROT are *inside* the cluster. In core
mode they are outside it, so the naive clamp would stop the window
dead at the element's own GAG.

Replace the clamp with an outward **walk**, governed by one principle
(O3, resolved):

> **A domain that survived `dante_filtering` is a real domain and
> blocks, unless it can be positively explained as part of this
> element.**

So the default verdict is *blocking*, and transparency has to be
earned. Starting at the seed and moving outward on each side, each
successive filtered domain is **transparent** only when **all** of the
following hold:

1. it is on the **same strand** as the seed;
2. it has LTR-retrotransposon support — `Final_Classification` starts
   with `Class_I|LTR`, or at least one `Region_Hits_Classifications`
   entry does (same test as §6.1);
3. its `Name` is listed in `accessory_5` / `accessory_3` for **this
   side of the core in element orientation** (so GAG 3' of a gypsy core
   blocks, and so does CHD 5' of it);
4. that accessory type has **not already been traversed** on this side
   — a second GAG going outward is the neighbouring element, not a
   second GAG of this one;
5. the accessory types encountered so far on this side appear in their
   canonical inward→outward order (gypsy 5': PROT then GAG; encountering
   GAG then PROT means the walk has crossed into another element);
6. where its own classification resolves to superfamily depth, that
   superfamily **agrees** with the seed's order-derived superfamily
   (O7, resolved). A domain that resolves no deeper than
   `Class_I|LTR` satisfies this vacuously — the rule only fires on a
   positive disagreement, never on missing information, so it costs
   nothing on the distant genomes where classification is shallow.

Anything else blocks: opposite strand, a non-LTR classification
(TPase, ENDO, HEL1/2 …), another seed's core domain, an accessory
domain on the wrong side or in the wrong order, one whose superfamily
contradicts the seed, and any domain that survived filtering but cannot
be placed by rules 1–6.

Note what this rule does **not** do: absence of accessory domains is
never blocking. A gypsy element with no detectable GAG or PROT simply
has nothing to traverse, and the walk runs to the offset cap. That is
the case core mode exists to handle, and it is unaffected by the
blocking-by-default polarity.

Rule 6 was initially left out, on the argument that GAG and PROT are the
least reliably classified domains on distant species and that requiring
agreement would reintroduce the classification dependency this mode
removes. Two things settle it the other way. First, the measurement:
**zero mismatches in 51 875 validated elements** across three genomes
(§5.2), so on everything we can check the rule never fires and costs no
sensitivity. Second, its asymmetry: because it triggers only on a
*positive* disagreement and treats an unresolved classification as
agreement, it cannot bite on the shallow-classification case that
motivated the objection. A GAG confidently called copia immediately 5'
of a gypsy core is evidence of an element boundary, and the walk should
stop there.

The count of such events is still reported per element
(`Accessory_Superfamily_Mismatch`) — now as a diagnostic of *why* a
window was truncated rather than a tolerated inconsistency — and feeds
§7's conflict reporting.

Stop at the first blocking neighbour, or at
`offset5prime`/`offset3prime` from the table, or at the sequence end,
whichever comes first. Reuse the existing `+100` idiom (extend 100 bp
into the blocking domain, since DANTE boundaries are approximate) and
the `offset2 = 300` overlap into the core, so the left and right windows
are built exactly as `get_ranges_left`/`get_ranges_right` build them,
only with per-seed limits substituted for `upstream_domain` /
`downstream_domain`.

Transparent domains traversed on the way out are recorded on the seed as
`Accessory_Domains` and reused in §7.

### 6.4 LTR, TSD and PBS

Unchanged: `get_TE(Lseq, Rseq, …, LTR_length = ltr_length)` →
`blast()` (tandem-repeat masking, `trim2TGAC`, TG…CA gate on both
copies, `-perc_identity 70`) → `evaluate_ltr()` → `get_best_ltr()`;
then `get_te_gff3()`, `add_pbs()`, `add_pbs_hemi()`, `get_te_rank()`.

`get_TE` already sorts BLAST hits by `qend − sstart` descending, i.e.
**innermost pair first** — precisely the "closest direct repeats"
requirement, and the reason a large window is tolerable.

Core mode adds structural gates on top of `get_best_ltr`'s existing
ones (TSD > 3 bp preferred, TE < 30 kb, LTR ≥ 100 bp), applied as a
post-filter so `ltr_utils.R` stays untouched:

- **G1** both LTRs lie entirely outside the core span (guaranteed by
  window construction, asserted anyway);
- **G2** `TE_Length ≥ core_span + 2 × min_ltr_length`;
- **G3** `TE_Length ≤ --max_te_length`;
- **G4** no **blocking** domain (§6.3 table) falls inside the called
  element — an element that swallows a Class II transposase or an
  opposite-strand domain is rejected, or demoted to `D`;
- **G5** the element contains exactly one seed. Two seeds inside one
  LTR pair means the outer repeat is not this element's LTR — reject
  the call outright (O4, resolved: reject, no nested-element reporting).
  Each seed is then re-examined on its own; the inner calls, if any,
  stand.
- **G6** with `--core_require_tsd`, reject `TSD == not_found`.

Elements failing a gate fall back to the rank-`D` track (§6.5) rather
than disappearing.

### 6.5 The rank-`D` track

Every filtered domain that does not end up inside a called element is
reported at rank `D`, whether or not it was ever core material (O5,
resolved). This keeps core-mode output comparable to lineage-mode
output and means nothing that survived filtering is silently dropped.

For v1 this reuses the lineage-mode machinery unchanged: build
`TE_partial` from the leftover domains with `get_domain_clusters_alt()`
+ `count_occurences_for_each_element()`, keep clusters with more than
one domain, then `trim_gr()` them against the called elements — exactly
the sequence already in `detect_putative_ltr.R`. No new code.

That machinery clusters on shared `Final_Classification`, which is the
very assumption core mode rejects, so on distant genomes the rank-`D`
grouping will be more fragmented than it should be. **Deferred**:
regrouping the leftover domains by domain *order* rather than
classification, consistent with §6.2. It does not affect core element
detection — the primary aim — and is tracked as follow-up work, not
part of this design.

---

## 7. Reannotation and classification

Only after boundaries are fixed. For each accepted element:

1. **Collect** every filtered DANTE domain inside the element on the
   element's strand — core, accessory, everything. These become the
   element's `protein_domain` children (as in lineage mode) and the
   evidence set below.
2. **Structural floor.** `Superfamily` from §6.2 is authoritative.
   Nothing below it may be reported without agreeing with it.
3. **Order-compatible lineages.** For each row of
   `lineage_domain_order.csv` under the called superfamily, test whether
   the observed domains appear in the lineage's canonical order — a
   subsequence test, with at most `--max_missing_domains` observed
   domains the lineage does not have at all. Missing domains are not
   penalised: an element carrying only RT/RH/INT genuinely *is*
   compatible with every gypsy lineage, and reporting that breadth is
   what `Lineage_Candidates` is for.

   This does **not** reuse `domain_distance()` (`ltr_utils.R:285`), as
   an earlier draft of this document proposed. That function computes
   `d_query_p == d_reference_p[d_reference_p %in% d_query_p]` without
   checking lengths, so whenever the query carries more domains than the
   reference the comparison recycles — R warns and the returned distance
   is meaningless. Lineage mode rarely hits this, because its clusters
   are already keyed to a single lineage; in core mode an incomplete or
   unexpected domain complement is the *normal* case, so the unsound
   path would be the common one. `domain_distance()` is left untouched.
   (Lineage mode can reach the same recycling when a cluster carries
   more domains than its lineage's reference order — worth a look, but
   out of scope here, since any change to it moves lineage-mode output.)
4. **Classification evidence.** For every collected domain, take
   `Final_Classification` plus every entry of
   `Region_Hits_Classifications`. Tally lineage-depth labels.
5. **Call.**
   - `Lineage_Call` is set only when a single lineage is both
     order-compatible and carried by a majority of the element's domains
     as their `Final_Classification`.
   - `Lineage_Candidates` is the support-sorted list of lineages that
     are order-compatible **and** appear in at least one domain's
     `Region_Hits_Classifications`.
   - `Lineage_Support` reports `n_supporting/n_total` domains.
6. **Demotion (`Final_Classification`).** Take the **lowest common
   ancestor** of the collected domains' `Final_Classification` values in
   the REXdb classification tree, then **clip it at the structural
   superfamily**. So:
   - all domains say `…/chromovirus/Tekay` → LCA is Tekay → report
     Tekay;
   - domains say Tekay and Chlamyvir → LCA is `…/chromovirus` → report
     chromovirus, with both in `Lineage_Candidates`;
   - domains disagree across superfamilies, or the classification-derived
     superfamily contradicts the order-derived one → report the
     order-derived superfamily and set `Classification_Conflict=true`.
     The order wins: it is structural evidence, the classification is a
     best-hit against a distant database.
   - `Classification_Demoted=true` whenever the reported depth is
     shallower than the deepest per-domain call.

This satisfies the requirement directly: *demote when necessary, report
the lineage when the evidence points to one, report the positive list
when ambiguous.*

---

## 8. Output

Feature types, ranks and file set are identical to lineage mode.
Additional attributes on the `transposable_element` feature:

```
source = dante_ltr_core
Final_Classification = Class_I|LTR|Ty3/gypsy        # §7.6
Name                 = Class_I|LTR|Ty3/gypsy        # kept == Final_Classification
Superfamily_Evidence = domain_order:RT,RH,INT
Core_Domains         = RT RH INT
Accessory_Domains    = PROT GAG CHD                 # element orientation
Accessory_Superfamily_Mismatch = 0                  # blocked the walk, §6.3
Lineage_Call         = Tekay                        # omitted if ambiguous
Lineage_Candidates   = Tekay,Chlamyvir,Reina
Lineage_Support      = 4/6
Classification_Demoted  = true
Classification_Conflict = false
Rank                 = DLTP
```

`Name`/`Final_Classification` carry the safe, possibly-superfamily-level
call, so `dante_ltr_summary`'s per-lineage tables and
`dante_ltr_to_library`'s clustering keep working; they simply see a
coarser label for elements where that is all the evidence supports.

`*_statistics.csv` keeps its exact shape (rows = classification,
columns = `D, DL, DLT, DLP, DLTP, RT_domain`) so the Python
`sum_up_stats_files()` merges chunk outputs unchanged.

Core-mode log additions (per chunk, aggregated by Python):

```
core candidates (RT/RH/INT)        : 5124   (412 via secondary LTR support)
ordered core seeds                 : 1387   (gypsy 1002 / copia 385)
seeds with LTR pair                : 1120
rejected by structural gates       :  118   (G2 31, G4 22, G5 65)
windows truncated by blocking rule :  944   (of 1387 seeds)
elements with lineage-level call   :  402
elements demoted to superfamily    :  718
classification conflicts           :   41
```

---

## 9. Risks

| risk | mitigation |
|---|---|
| Large windows → spurious direct repeats | innermost-pair ordering, TG…CA on both copies, TSD preference, tandem masking, G1–G3; measured offsets (§5.2) are ~half the first estimate, so the exposure is smaller than feared |
| Two adjacent elements merged into one | G5 (one seed per element); blocking-by-default walk stops at the neighbour's core, at a repeated accessory type, or at an out-of-order one |
| Non-LTR RT (LINE, pararetrovirus) seeding false elements | RT alone never seeds; the ordered triplet with INT is required, and LINEs have no INT |
| Superfamily mis-call from a mis-ordered fragment | order requires all three domains on one strand within `core_max_span`; a fragment yields no seed, not a wrong seed. Measured 0/1135 discordant (§5.2) |
| Core mode quietly worse than lineage mode on good genomes | it is opt-in; §10 requires a concordance test on REXdb-covered data |
| **Sensitivity capped by the domain filter, not the algorithm** | up to ~24 % of validated elements lack a complete filtered core (§5.2, Alyr; 2–7 % elsewhere). Partly mitigable via `--min_relative_length_core` (§6.0.2), worth ~4–10 % more seeds (§5.3, §5.4). The residual bound must be stated in the README so users do not read core mode as strictly more sensitive on well-covered genomes |
| Blocking-by-default truncates windows on domain-dense regions | §6.3's five transparency rules cover the legitimate cases; §10 tests each. Elements lost this way still appear at rank `D` (§6.5) |
| `Region_Hits_Classifications` widening the candidate set too far | it is only used for `Lineage_Candidates` (advisory) and for secondary LTR support, never for `Final_Classification` |

---

## 10. Validation and tests

**Concordance (specificity).** Run both modes on a REXdb-covered
genome. The target must be **conditional on the core being present**,
because §5.2 measured that 2–25 % of validated elements (genome
dependent; Alyr is the worst) do not carry all three core domains
through `dante_filtering` — core mode cannot find
the rest, by construction, and an unconditional target would fail for a
reason that has nothing to do with the algorithm:

> Of the lineage-mode `DLT`/`DLTP` elements **whose three core domains
> survive filtering**, core mode must recover ≥ 90 % with boundaries
> within ±20 bp.

Report both numbers — conditional recovery and the unconditional
fraction — so the filter-imposed ceiling stays visible and is not
mistaken for an algorithmic loss. The excess (core-mode calls with no
lineage-mode counterpart) is the specificity signal: a large excess on
a well-covered genome means the structural gates of §6.4 are too loose.

**Superfamily correctness.** For every element found by both modes,
compare the order-derived superfamily against lineage mode's
classification. §5.2 measured 1135/1135 agreement on
`test_data/g1_dante_ltr.gff3`, so this should be ≈0 discordant; treat
any regression as a bug in seeding, not a tolerance to widen.

**Sensitivity.** Run both modes on Drapa, where lineage mode currently
returns **zero** complete elements against 1 984 ordered core seeds
(§5.4). The seed count is an upper bound; the test measures how many
become elements with an LTR pair, and checks their length and
LTR-identity distributions for plausibility. Compare against
`--fallback_mode coarse2/coarse3` on the same input, since that is the
alternative a user would otherwise reach for (§11) — and which, on a
genome with no complete domain complements, should also find little.

**Filter sensitivity.** Sweep `--min_relative_length_core` end-to-end on
Drapa and Pisum, reporting *elements* (not just seeds) at each step.
§5.3/§5.4 measured the seed response; what remains unmeasured is whether
the extra seeds convert to elements with LTRs or are simply discarded at
§6.4. That is the number that justifies the 0.3 default.

**Unit-ish tests** (run via a new `tests/core.sh` wired into `tests.sh`,
PR tier):
- seeding on a synthetic GFF3: tandem array → two seeds; minus strand;
  RT+RH only → no seed; `core_max_gap` exceeded → no seed; duplicated
  RT within one run.
- blocking rule (§6.3): opposite-strand domain, TPase, neighbour core,
  GAG 3' of a gypsy core, second GAG going outward, GAG-before-PROT
  order, copia-classified GAG 5' of a gypsy core → all truncate the
  window; PROT then GAG → window passes through; a GAG resolving only to
  `Class_I|LTR` → passes through (rule 6 must not fire on missing
  information); no accessory domains at all → window runs to the offset
  cap.
- G1–G6 gates, one synthetic case each.
- LCA demotion (§7.6): the four cases.
- `--mode core` end-to-end smoke on `sample_genome_part.fasta`.
- fragment merging (§6.0.1): synthetic split RT with complementary
  `Best_Hit_DB_Pos` merges; two genuine copies with overlapping
  `Best_Hit_DB_Pos` do not; a split beyond `--split_max_gap` does not;
  merging happens before filtering, so a pair of sub-threshold fragments
  is rescued.
- `calibrate_core_constraints.py` on `test_data/g1_dante_ltr.gff3`
  reproduces the §5.2 numbers — this pins the measurement that the
  shipped table is derived from.

---

## 11. Relationship to the other approaches

| | requires full domain complement | requires lineage-depth classification | detects superfamily from |
|---|---|---|---|
| lineage mode (default) | yes | yes | classification |
| `--fallback_mode coarse2/3` | yes | no (demoted) | classification |
| `--mode core` | **no** (RT/RH/INT only) | **no** | **domain order** |

Core mode is the strictly more permissive structural criterion and the
only one that survives a genome where GAG and PROT are undetectable.
All three stay available.

---

## 12. Resolved decisions

| | question | resolution |
|---|---|---|
| **O1** | lower `min_similarity` for core mode? | **No.** `Similarity` is *exactly* non-binding -- identical seed counts at 0.20/0.30/0.40 on both Pisum and Drapa (§5.3, §5.4). Keep 0.4 and expose it. The real lever is `Relat_Length`, relaxed for core domains only via `--min_relative_length_core` (default 0.3). §3.1, §6.0.2 |
| **O8** | re-measure fragmentation and thresholds on an uncovered genome | **Done** (§5.4, Drapa). Both expectations refuted: the frameshift excess vanishes into background, and the `Relat_Length` gain is *smaller* than on Pisum. Defaults `--min_relative_length_core 0.3` and `--split_max_gap 50` stand. |
| **O3** | ambiguous domain in the window: transparent or blocking? | **Blocking.** Surviving the filter means the domain is real; transparency must be positively earned. §6.3 |
| **O4** | two seeds in one LTR pair: reject or report as nested? | **Reject.** No nested-element reporting. §6.4 G5 |
| **O5** | rank-`D` track from non-core clusters too? | **Yes** — anything that passes the filter and is not inside an element is reported at rank `D`. Order-aware regrouping of that track is deferred; core element detection is the primary aim. §6.5 |
| **O6** | does core-mode output feed `dante_ltr_solo`? | **Follow-up**, out of scope here. |
| **O2** | seed scoring function | **Shortest span, ties on `sum(Similarity)`.** No free parameter; the measured span distribution (q99 ≈ 2.9–3.6 kb across four genomes) supports the compactness prior. §6.2 step 3 |
| **O7** | must a traversed accessory domain's superfamily agree with the seed's? | **Yes -- disagreement blocks** (§6.3 rule 6). Fires only on positive disagreement, so a shallowly-classified domain passes vacuously; measured 0 mismatches in 51 875 elements, so it costs no sensitivity where checkable. |

## 13. Nothing open

Every question raised during design review is resolved (§12). What
remains before implementation is measurement, not decision:

- **The 0.3 default for `--min_relative_length_core` is justified on
  seed counts, not element counts.** §5.3/§5.4 measured how many extra
  *seeds* the relaxation buys (+4–10 %); whether those seeds survive the
  LTR search of §6.4 and become elements is untested, because it needs a
  working implementation. §10 names this as the test that settles the
  default, and it should be run before the value is frozen in the
  shipped table.
- **The constraints table is built from four genomes** (§5.2, §5.4).
  The Darwin runs should widen it — `utils/calibrate_core_constraints.py`
  per genome, taking the maximum of the per-genome offsets.

## 14. Follow-up work (explicitly out of scope)

- Regrouping the rank-`D` track by domain order instead of shared
  classification (§6.5).
- Feeding core-mode elements to `dante_ltr_solo` (O6).
- Any change to `dante_ltr_to_library`, `dante_ltr_summary` or
  `clean_ltr.R` — they consume core-mode output unchanged (§3.2).
