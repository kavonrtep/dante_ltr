# Implementation plan — CI tests + automated conda release

## 1. Goal

Replace the current manual workflow
(`kavonrtep/recipes` → local Singularity → `anaconda login` → `conda mambabuild`)
with tag-driven GitHub Actions:

1. **Every push / PR** runs a smoke + short test suite (~3 min total).
2. **A release tag push** runs the long test as a gate, then builds and
   uploads the conda package to the `petrnovak` anaconda channel.
3. No manual recipe edits, no manual sha256, no interactive login.

Own channel only for now; bioconda submission is a later task once the
release flow is stable.

## 2. Test layout

```
tests/
├── data/
│   ├── smoke/
│   │   ├── genome.fasta      # 1 contig, ~kb–10 kb, one complete TE region
│   │   └── dante.gff3        # DANTE lines for that contig only
│   └── short/
│       ├── genome.fasta      # ~5 contigs, several lineages, DLTP+DLT+DLP
│       └── dante.gff3        # DANTE lines for those contigs only
├── smoke.sh                  # < 30 s — --help/--version, tiny pipeline, asserts
├── short.sh                  # 1–3 min — full pipeline, asserts non-trivial counts
└── long.sh                   # 10–30 min — on existing test_data/g1.*
```

Root `tests.sh` becomes a thin dispatcher:

```bash
./tests.sh smoke   # default
./tests.sh short
./tests.sh long
./tests.sh all
```

Backwards-compat: when called with a number, treat as CPU count for
`long.sh` (matches the current behaviour of `./tests.sh 4`).

### Subset construction

Full DANTE data live in `test_data/sample_DANTE.gff3` +
`test_data/sample_genome.fasta`; `tmp/test_output1.gff3` is a previous
`dante_ltr` run whose `transposable_element` features tell us which
contigs are interesting.

**Principle**: *contig-level slicing only*. We pick whole contigs; no
coordinate remapping. This keeps the DANTE GFF3 coordinates valid
without transformation and avoids partial-element / cross-contig
pathologies.

Selection algorithm (run once, captured in a short R script so it's
reproducible):

1. Load `tmp/test_output1.gff3`; keep only `transposable_element` rows.
2. Group by `seqnames`; compute per-contig rank mix and lineage mix.
3. **smoke contig**: smallest contig that contains ≥ 1 DLTP element and
   whose length is ≤ 100 kb. One contig only.
4. **short contigs**: greedy — add contigs until we have at least 3
   distinct lineages, at least 1 DLTP and 1 DLT, at least 1 solo LTR
   will likely be detectable (i.e. multiple copies of ≥ 1 lineage), and
   total length < 5 Mb.
5. For each selected contig:
   - slice the FASTA record out
   - `grep -P "^<contig>\t"` from `sample_DANTE.gff3` into the subset GFF3
   - (optional) prepend the GFF3 header lines

Output goes to `tests/data/smoke/` and `tests/data/short/`. Expected
sizes: smoke < 100 kb FASTA, short < 5 MB FASTA. Both commit fine.

### Assertions

Each `*.sh` produces output in `tmp/tests/<level>/` and validates:

**smoke.sh**
- every executable prints `--help` and `--version` without error
- `dante_ltr` on the subset exits 0 and produces a non-empty GFF3
- produced GFF3 parses as valid GFF3 (first line `##gff-version 3`)
- at least 1 `transposable_element` feature

**short.sh**
- all smoke checks
- ≥ 1 `transposable_element` of `Rank=DLTP` or `DLT`
- `dante_ltr_to_library` exits 0 and produces non-empty FASTA
- `dante_ltr_solo` exits 0; `solo_ltr.gff3` non-empty;
  `solo_ltr_statistics.csv` has ≥ 1 lineage row
- `solo_ltr_raw.gff3` has ≥ `solo_ltr.gff3` feature count (rep
  collapsing works)

**long.sh**
- runs the existing `test_data/g1*` end-to-end
- compares lineage counts against a loose baseline (ranges, not exact
  values, so minor version-to-version drift doesn't fail CI)

All three use `set -euo pipefail`; on any assertion failure exit with
non-zero and print the expected vs observed values.

## 3. Conda recipe in-repo

```
conda/
└── dante_ltr/
    ├── meta.yaml
    └── build.sh
```

**meta.yaml** — ported from `kavonrtep/recipes/dante_ltr/meta.yaml`:

- `{% set version = load_file_data('../../version.py')['__version__'] %}`
  → read directly from `version.py` so the tag and the package always
  agree (also catches the "forgot to bump version.py" case).
- `source: path: ../..` — **local path**, no URL, no sha256. CI checks
  out the repo at the tag, the build uses that checkout.
- `requirements.run:` full list from `requirements.txt` (adds
  `r-dplyr=1.0.7` ... and **`mafft>=7.490`** which is missing from the
  current recipe even though it's in `requirements.txt`).
- `test.commands:` `dante_ltr --help`, `dante_ltr --version`,
  `dante_ltr_solo --help`, `dante_ltr_to_library --help`.
- `about` / `extra` unchanged.

**build.sh** — same pattern as today:

```sh
#!/bin/sh
set -x -e
DANTE_LTR_DIR=${PREFIX}/share/dante_ltr
mkdir -p ${PREFIX}/bin ${DANTE_LTR_DIR}
cp -r . ${DANTE_LTR_DIR}
ln -s ${DANTE_LTR_DIR}/dante_ltr            ${PREFIX}/bin/dante_ltr
ln -s ${DANTE_LTR_DIR}/dante_ltr_solo       ${PREFIX}/bin/dante_ltr_solo
ln -s ${DANTE_LTR_DIR}/dante_ltr_to_library ${PREFIX}/bin/dante_ltr_to_library
ln -s ${DANTE_LTR_DIR}/dante_ltr_summary    ${PREFIX}/bin/dante_ltr_summary
ln -s ${DANTE_LTR_DIR}/clean_ltr.R          ${PREFIX}/bin/clean_ltr.R
```

Adds `dante_ltr_solo` (missing from the current recipe). Uses `cp -r .`
instead of `cp -r *` to also copy dotfiles if any.

## 4. GitHub Actions workflows

### `.github/workflows/tests.yml` — on push and PR

```yaml
name: tests
on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  smoke-and-short:
    runs-on: ubuntu-latest
    defaults: { run: { shell: bash -el {0} } }
    steps:
      - uses: actions/checkout@v4

      - uses: conda-incubator/setup-miniconda@v3
        with:
          miniforge-version: latest
          activate-environment: dante_ltr
          channels: conda-forge,bioconda,petrnovak
          channel-priority: strict

      - name: Install runtime dependencies
        run: mamba install -y --file requirements.txt

      - name: Smoke test
        run: ./tests.sh smoke

      - name: Short test
        run: ./tests.sh short
```

Key points:
- `miniforge-version: latest` matches the Singularity base (mambaforge).
- Channel priority + order: conda-forge, bioconda, petrnovak (same as
  local Singularity).
- Install straight from `requirements.txt` — same file the manual
  install uses.
- Smoke and short in one job (caching dominates — splitting would
  double the ~2 min env-setup cost).

### `.github/workflows/conda-release.yml` — on tag push

```yaml
name: conda-release
on:
  push:
    tags: ['[0-9]+.[0-9]+.[0-9]+*']

jobs:
  build:
    runs-on: ubuntu-latest
    defaults: { run: { shell: bash -el {0} } }
    steps:
      - uses: actions/checkout@v4

      - name: Assert tag matches version.py
        run: |
          V=$(python -c "exec(open('version.py').read()); print(__version__)")
          [ "$V" = "$GITHUB_REF_NAME" ] || {
            echo "version.py ($V) does not match tag ($GITHUB_REF_NAME)"; exit 1; }

      - uses: conda-incubator/setup-miniconda@v3
        with:
          miniforge-version: latest
          channels: conda-forge,bioconda,petrnovak
          channel-priority: strict

      - name: Install build tools + runtime deps
        run: |
          mamba install -y conda-build boa anaconda-client
          mamba env create -n test_env --file requirements.txt

      - name: Long test (release gate)
        run: |
          conda activate test_env
          ./tests.sh long

      - name: Build conda package
        run: |
          conda mambabuild -c conda-forge -c bioconda -c petrnovak \
            --output-folder build_out conda/dante_ltr

      - name: Upload to anaconda
        env:
          ANACONDA_API_TOKEN: ${{ secrets.ANACONDA_API_TOKEN }}
        run: |
          anaconda -t "$ANACONDA_API_TOKEN" upload \
            --user petrnovak \
            --label main \
            build_out/noarch/dante_ltr-*.tar.bz2
```

Gate structure:
1. tag-name ↔ `version.py` consistency check (fail fast if forgotten).
2. long test must pass — guarantees a broken commit never reaches
   anaconda.
3. build.
4. upload with token; user fails gracefully if secret isn't set.

## 5. Secrets

One repo secret: `ANACONDA_API_TOKEN`.

Generation at `anaconda.org` → *Settings* → *Access* → *Create token*:
- write scope: ✔ Allow write access to the API site
- no read scope needed
- 1-year expiry; add a reminder to rotate

Added at `github.com/kavonrtep/dante_ltr` → *Settings* → *Secrets and
variables* → *Actions* → *New repository secret*. Name exactly
`ANACONDA_API_TOKEN`.

## 6. Migration notes

- `recipes/dante_ltr/` in `kavonrtep/recipes` is not touched by this
  change. Once the GHA flow publishes successfully a couple of times I'd
  replace that folder's contents with a short README pointing to the
  in-repo recipe; keep the other packages (`dante`, `dante_tir`, …) as
  today.
- `Singularity.def` stays where it is — still useful for manual debug
  builds on your HPC.
- Existing `tests.sh` at the repo root becomes a dispatcher; the old
  invocation `./tests.sh 4` maps to `./tests.sh long 4` for backwards
  compatibility.

## 7. Local verification before pushing

Two checks must pass before we enable CI:

1. `./tests.sh smoke` locally → exit 0.
2. `./tests.sh short` locally → exit 0.
3. A local `conda mambabuild conda/dante_ltr` in a fresh env → succeeds
   and passes the recipe `test.commands`. (Optional but recommended.)

If any of these fail we refine the subset data / assertions before
committing.

## 8. Explicitly out of scope (deferred)

- Submitting to `bioconda-recipes`.
- Multi-platform conda builds (osx-arm64, linux-aarch64 — currently
  noarch is fine since everything is scripts + R/python via conda).
- Caching the conda env between workflow runs (can halve CI time once
  everything is green — optimisation later).
- Publishing to `Anaconda.org` under a second label for beta releases.
