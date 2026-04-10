# Solo LTR Detection — Design Document

## Biological Background

### What is a solo LTR?

LTR retrotransposons (LTR-RTs) are flanked on both ends by identical Long Terminal Repeats (LTRs). Over evolutionary time, the two LTRs of a single element can undergo intrachromosomal homologous recombination, looping out the internal coding sequence and leaving behind a single "solo" LTR. Solo LTRs therefore:

- Lack all protein-coding domains (no RT, INT, GAG, etc.)
- Retain the LTR sequence (and structural signals: TG...CA termini, PBS-adjacent position on the 3' side of the 5' LTR)
- Are flanked by **Target Site Duplications (TSDs)** inherited from the original insertion event

In many plant genomes, solo LTRs outnumber complete elements by 5–10-fold. They represent a major component of repetitive DNA that is missed by domain-based annotation tools like DANTE.

### Why TSD is the key distinguishing feature

Any region of genomic DNA that shares sequence similarity to an LTR could arise from:
1. A *bona fide* solo LTR (genuine insertion remnant)
2. An internal fragment of a nested or fragmented complete element not captured by DANTE_LTR
3. A spurious low-complexity match

Only scenarios 1 and 3 produce a hit that sits entirely outside any annotated complete element, but scenario 3 will lack TSDs. The TSD acts as a molecular fossil of the original integration: at insertion, the integrase creates a staggered cut and fills in both ends with identical short sequences. After recombination, both TSDs survive flanking the solo LTR. A confirmed TSD pair thus unambiguously marks a real insertion event.

Scenario 2 (unannotated fragment) is subtler. Such a fragment will typically carry conserved sequence extending beyond the LTR boundary — a 5'UTR or PBS signature immediately 3' of a 5'-LTR-like hit, or internal sequence immediately 5' of a 3'-LTR-like hit — and it will lack TSDs.

### TSD length estimation

TSD lengths are lineage-specific (typically 4–6 bp) and are **estimated empirically** from the input DANTE_LTR annotation. DLT and DLTP ranked elements already carry validated TSD sequences in the GFF3. By grouping these per lineage and taking the modal TSD length, we obtain a reliable per-lineage prior that is more robust than a hardcoded table, because not all elements will have a detected TSD and occasional errors exist in individual calls. Where insufficient DLT/DLTP data exist for a lineage, fall back to scanning 4–6 bp.

---

## Prior Work — NGS-based solo LTR Detection

A related approach was previously used to detect and quantify solo LTRs from NGS read data (Supplementary Figure S4, Novak et al. 2015, PLoS ONE). Although that method operated on short reads rather than an assembled genome, its core logic informs the present design.

**Key principle**: The junction between the LTR 3' end (ending in `...CA`) and the 5'UTR is the structural hallmark of a complete element. By extracting 30 bp sequence tags from both the **LTR_3'end** and the **5'UTR** region of known complete elements, two tag databases are built. Reads (or, in our context, genomic windows) that hit the LTR_3'end tag but *lack* a hit to the 5'UTR tag represent LTR/insertion-site junctions — i.e., solo LTRs or the 3' LTR of a complete element. Reads hitting *both* tags are LTR/5'UTR junctions belonging to complete elements.

**Rsf ratio**: The ratio of solo LTRs to complete elements per lineage is expressed as:

```
Rsf = (Lx − LU) / LU
```

where *LU* = number of LTR/5'UTR junction reads and *Lx* = number of LTR_3'end-only hits (insertion-site reads). In the genome annotation context this translates to a simple count ratio per lineage reported in the statistics output.

**Adaptation for assembled genome annotation**: Rather than 30 bp tags on short reads, we work with full LTR sequences on an assembled genome, but the 5'UTR junction check is directly reusable as a discriminator step (see Step 4b below).

---

## Problem Formulation

**Input:**
- DANTE_LTR GFF3 annotation of a genome (`long_terminal_repeat`, `target_site_duplication`, and `transposable_element` features; ranks DL / DLT / DLP / DLTP). DLT and DLTP elements additionally provide empirical TSD lengths per lineage.
- Reference genome FASTA

**Goal:**
- Identify all solo LTR insertions genome-wide that are not already part of an annotated complete element
- Report each candidate with: genomic coordinates, classification (lineage), % identity to the reference LTR, TSD sequence and coordinates
- Separately flag hits that lack TSDs and carry signatures of incomplete/fragmented elements
- Report the Rsf ratio per lineage as a summary statistic

**Acceptance criteria:**

| Rank         | Criteria |
|--------------|----------|
| `SL`         | Similarity hit outside complete elements + confirmed TSD pair |
| `SL_noTSD`   | Similarity hit outside complete elements, no TSD — reported separately, not classified as solo LTR |

---

## Algorithmic Design

### Step 1 — Build a non-redundant LTR reference library

Extract all `long_terminal_repeat` child features from the DANTE_LTR GFF3. Prefer higher-rank elements (DLTP > DLT > DLP > DL) as their LTR boundaries are more precisely determined (TSD validation in DLTP/DLT provides independent confirmation of boundary accuracy).

For each extracted LTR, also retrieve a short flanking window (15 bp on each side). The extended sequences are used only for the MSA step; the pure LTR sequences are used as BLAST queries.

**Consensus construction with MAFFT**: For each cluster, compute a multiple sequence alignment with MAFFT using the LTR sequences including flanking. The MSA reveals:
- The true TG...CA boundary positions even when individual element boundaries are off by 1–2 bp
- The TSD positions in the flanking context (often visible as a conserved short block immediately outside the TG/CA)
- The 5'UTR tag region immediately 3' of the LTR (see Step 4b)

The MSA consensus becomes the primary BLAST query for that cluster, replacing any single representative. This improves sensitivity for diverged solo LTRs.

Cluster with MMseqs2 at **90 % identity** (reusing `utils/mmseq_clustering.R`). Per-lineage clustering is preferred over genome-wide clustering to avoid merging related lineages.

**5'UTR tag database (auxiliary)**: From each complete element in the GFF3, extract the 30 bp immediately downstream (3') of the 5' LTR end. These tags encode the LTR/5'UTR junction and are used in Step 4b to discriminate solo LTRs from unannotated complete elements. Store tags with lineage labels, analogous to the NGS approach in the prior work.

Output:
- `<prefix>_LTR_library.fasta` — consensus LTR sequences per cluster, headers encoding lineage
- `<prefix>_5UTR_tags.fasta` — 30 bp 5'UTR junction tags per lineage for Step 4b

### Step 2 — Genome-wide similarity search

Run BLAST (`blastn`) of the LTR library against the genome:

```
blastn -task blastn -query LTR_library.fasta -db genome.fasta \
  -perc_identity 80 -dust no \
  -outfmt "6 qaccver saccver pident length qlen qstart qend sstart send sstrand evalue bitscore"
```

Key filter parameters applied to raw output:
- Minimum identity: 80 %
- Minimum alignment length: 80 % of query LTR length (`length / qlen ≥ 0.8`)
- Maximum e-value: 1e-5

BLAST is preferred over MMseqs2 for its sensitivity and directly compatible output format (same format used throughout the existing pipeline). For very large genomes, the BLAST database is built once from the full genome FASTA; the library is small. The chunk-splitting and CPU-distribution strategy from `dante_ltr` (Python wrapper) is reused for parallelisation.

### Step 3 — Remove hits overlapping known complete elements

Intersect BLAST hits with all features from the input GFF3 (including child features to catch LTR sequences internal to annotated elements) using `bedtools intersect` (reusing `trim_gr()` from `utils/ltr_utils.R`). Discard any hit overlapping an annotated feature by more than a small tolerance (20 bp).

Surviving hits include both genuine solo LTRs and fragments of unannotated complete elements; Steps 4a and 4b discriminate between them.

### Step 4 — TSD validation and fragment discrimination

#### Step 4a — TSD check

For each surviving hit, mirrors the logic in `evaluate_ltr()` in `utils/ltr_utils.R`:

```
genome: ...[ TSD_L ][ ←  solo LTR hit  → ][ TSD_R ]...
              4–6 bp                          4–6 bp
```

1. Determine expected TSD length for the hit's lineage from the empirical map derived from DLT/DLTP annotations (fall back to scanning 4–6 bp if data are insufficient).
2. Extract expected-length bp immediately 5' of the hit start → `TSD_L_seq`.
3. Extract expected-length bp immediately 3' of the hit end → `TSD_R_seq`.
4. Compare (same logic as `evaluate_ltr()`):
   - Exact match → TSD confirmed
   - ≤1 mismatch and length ≥ 5 bp → TSD confirmed with mismatch
   - No match → no TSD; proceed to Step 4b
5. Minus-strand hits: reverse-complement both sequences before comparison.

#### Step 4b — 5'UTR junction check and PBS check (no-TSD hits only)

Adapted from the NGS-based prior work. For each hit that fails the TSD check:

Two junction checks cover both LTR boundaries, plus a PBS check as secondary confirmation:

**5'UTR junction check** (3' boundary): Extract 30 bp immediately 3' of the hit end (strand-aware). BLAST against the `5UTR_tags` database. A match means the hit is the **5' LTR** of an unannotated element still attached to its 5'UTR/PBS — the same LTR_3'end/5'UTR junction signature used in the NGS-based prior approach.

**PPT junction check** (5' boundary): Extract 30 bp immediately 5' of the hit start (strand-aware). BLAST against the `PPT_tags` database (30 bp immediately upstream of each annotated 3' LTR, capturing the PPT and flanking 3'UTR). A match means the hit is the **3' LTR** of an unannotated element with internal sequence still present upstream.

**PBS check**: BLAST the same 30 bp 3' window against the tRNA database (reusing `add_pbs()` / `add_pbs_hemi()`). Secondary confirmation for the 5'-LTR scenario.

All three are diagnostic attributes on `SL_noTSD` features, not hard filters. A hit negative for both junction checks is more likely a genuine old insertion where the TSD has diverged beyond recognition.

### Step 5 — Output

**GFF3** — one feature per hit, with `target_site_duplication` child features for confirmed solo LTRs:

```
##gff-version 3
Chr1  dante_ltr  solo_LTR  1000  1550  .  +  .  ID=soloLTR_00000001;Final_Classification=Class_I|LTR|Ty3/gypsy|chromovirus|CRM;Identity=91.3;Coverage=0.95;TSD=TTGCA;Rank=SL
Chr1  dante_ltr  target_site_duplication  995   999   .  +  .  Parent=soloLTR_00000001
Chr1  dante_ltr  target_site_duplication  1551  1555  .  +  .  Parent=soloLTR_00000001

Chr1  dante_ltr  solo_LTR  3000  3480  .  -  .  ID=soloLTR_00000002;Final_Classification=Class_I|LTR|Ty1/copia|Ivana;Identity=83.1;Coverage=0.82;TSD=not_found;Rank=SL_noTSD;UTR5_junction=negative;PBS_check=negative
```

Attribute fields:
- `Final_Classification` — lineage of the best-matching library LTR
- `Identity` — % nucleotide identity to the best-matching library LTR
- `Coverage` — fraction of library LTR length covered by alignment
- `TSD` — TSD sequence, or `not_found`
- `Rank` — `SL` (confirmed solo LTR) or `SL_noTSD` (unconfirmed)
- `UTR5_junction` — `positive` / `negative` / `not_applicable` (SL_noTSD only; positive = hit is likely the 5' LTR of an unannotated element)
- `PPT_junction` — `positive` / `negative` / `not_applicable` (SL_noTSD only; positive = hit is likely the 3' LTR of an unannotated element)
- `PBS_check` — `positive` / `negative` / `not_applicable` (SL_noTSD only; secondary confirmation for 5'-LTR scenario)

**Statistics CSV** — counts of SL and SL_noTSD per lineage plus the Rsf ratio (SL / complete elements from DANTE_LTR input), analogous to existing `_statistics.csv`.

---

## Integration with Existing Pipeline

New standalone script `dante_ltr_solo`, callable after `dante_ltr`:

```
dante_ltr_solo -g DANTE_LTR_annotation.gff3 \
               -s genome.fasta \
               -o solo_output \
               [-c CPU] \
               [-i min_identity 80] \
               [-C min_coverage 0.8]
```

### Reuse of existing components

| Existing component | Reuse |
|--------------------|-------|
| `utils/extract_fasta.R` | Extract LTR sequences and flanking by coordinates |
| `utils/mmseq_clustering.R` | Cluster LTRs at 90 % identity per lineage |
| `utils/ltr_utils.R::trim_gr()` | Remove hits overlapping annotated elements |
| `evaluate_ltr()` TSD logic | Port directly for TSD validation (Step 4a) |
| `add_pbs()` / `add_pbs_hemi()` | PBS check for SL_noTSD hits (Step 4b) |
| `dante_ltr` chunk-splitting logic | Parallelise genome-wide BLAST search |
| `Gff3Feature` class in `dante_ltr` | GFF3 parsing in the Python wrapper |

---

## Open Questions / Design Decisions

1. **MSA consensus library**: MAFFT will be used for consensus construction in Step 1. This is included in the initial implementation. MAFFT needs to be added as a conda dependency.

2. **SL_noTSD reporting**: These hits cannot be classified as solo LTRs. They are always reported (with 5'UTR junction and PBS diagnostic attributes) to allow users to distinguish fragmented complete elements from edge cases where TSD has genuinely diverged beyond recognition.

3. **Nested elements**: A solo LTR inside another unannotated TE will survive Step 3. This is a very edge case depending on DANTE_LTR completeness and is deferred to future work.

4. **BLAST for search**: BLAST is preferred over MMseqs2 for sensitivity and interpretability. It is already a dependency and shares output format with the PBS-check BLAST calls. MMseqs2 can be added as an optional fast backend in a future version if runtime becomes a bottleneck on very large genomes.
