# **DArTseq_InDel_Finder**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23775495.svg)](https://doi.org/10.5281/zenodo.23775495)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---
### InDels that Differentiate Samples (Cultivars / Species)

A pipeline to detect, filter, and select **InDels** obtained from DArTseq that
differentiate groups of samples (e.g., cultivars or species). Beyond detecting the
InDels, it designs HRM primers and predicts, per genotype, the High Resolution Melting
separation — so the output is a set of prioritized, ready-to-test HRM assays.

> It is the evolution (InDel + HRM) of [`DArTseq_SNP_Finder`](https://github.com/BioinfoUfes/DArTseq_SNP_Finder),
> developed by the same group (BioInfo/UFES).

From the DArTseq reports and one (or more) reference genome(s) of your species, the tool:

1. **infers indels** by comparing tags that fall at the same genome locus (DArTseq-SNP
   only reports substitutions; indels are inferred);
2. **maps and checks the specificity** of each tag in the genome you provide (works with
   any assembly, including separate haplotypes);
3. **designs HRM primers** (Primer3) with a short amplicon flanking the indel;
4. **predicts each sample's melting curve from its genotype** (A = long homoduplex,
   B = short homoduplex, AB = heterozygote with heteroduplex), via a local model
   (Zimm-Bragg/Ising with nearest-neighbor parameters) — without relying on any website;
5. **builds decision trees** that identify the samples in the fewest analyses, prioritizing
   real HRM separation (difference-plot peak).

It is **generic**: it works for any species — just swap the reports and the genome.

### What's new in v2.0.0

- Proper **affine-gap DNA aligner** (Biopython) for InDel calling, with a difflib fallback.
- **Exact-match in-silico PCR** screen of each primer pair against the reference (products
  counted per haplotype).
- **Honest `N` handling**: a double-null (`N`) is never treated as allele `B`; pairs that
  rely on `N` are flagged, never reported as STRONG.
- **Honest HRM labeling**: difference-plot peaks are uncalibrated, relative predictions —
  not instrument-calibrated values or probabilities of success.
- **Machine-readable outputs** (TSV + JSON) next to the XLSX, plus unit tests and CI.

See `CHANGELOG.md` for the full list.

---

## Installation

With conda (recommended):

```bash
conda env create -f environment.yml
conda activate indel-hrm-finder
```

Or with pip:

```bash
pip install -r requirements.txt
```

Dependencies: `primer3-py`, `biopython`, `numpy`, `openpyxl`, `matplotlib`.
Neither BLAST nor anything else needs to be installed.

## Inputs

- **DArTseq reports** in a folder (e.g. `reports/`): one `*SNP_mapping*.csv` **and** one
  `*SilicoDArT*.csv` from your DArT order.
- **Reference genome(s)** of your species, in FASTA (`.fa`/`.fa.gz`). You may provide
  more than one file (e.g. two haplotypes of a phased assembly). Localization is done by
  the tag sequence itself, so **any assembly works**.
- *(optional)* **Name CSV** for the samples, with columns `code,name`, so the reports
  come out with readable names instead of the codes (see `examples/`).

## Usage

### Interactive mode (asks the questions)

```bash
python run.py
```

It asks for: the reports folder, the genome(s), the name CSV (optional), the mode
(tree or panel), which samples to compare, your kit's HRM conditions (Mg2+ and Mono+)
and the output folder.

### Command-line mode

```bash
python run.py \
  --reports reports \
  --genome genomeHapA.fa.gz genomeHapB.fa.gz \
  --metadata examples/example_samples.csv \
  --samples 7011,7041,7115,7167 \
  --mode tree \
  --melt-mg 3 --melt-mono 20 \
  --out output
```

- `--mode tree` builds the **decision tree** for the group of samples given in
  `--samples` (omit to use all).
- `--mode panel` builds the **minimal panel** that distinguishes all the samples.
- `--melt-mg` / `--melt-mono`: use the concentrations of **your** HRM master mix.

## Outputs (in the `--out` folder)

- `HRM_markers.xlsx`
  - **Separation_per_pair** — for each pair of samples, the best markers and the
    predicted difference-plot peak (your tool for choosing primers).
  - **Trees** — the decision trees (best first), with the step-by-step flow.
  - **Primers** — primers (F/R, Tm), amplicons, mapping hits and in-silico PCR products.
  - **Method** — parameters and notes.
- `candidates.tsv`, `primers.tsv`, `per_pair.tsv` — the same information as plain tables,
  for scripting and sharing.
- `run_summary.json` — parameters, counts and the chosen tree(s), machine-readable.
- `curves/` — one image per marker: the predicted curve of **each sample by genotype** +
  the difference plot.
- `amplicons.fasta` — long/short alleles of each marker, for checking in uMelt.

## Testing

```bash
pip install pytest
pytest -q
```

The tests use small synthetic fixtures and run in seconds. They verify software behavior
(alignment, genotype states, `N` handling, in-silico PCR, panel selection), not biological
validity on real data.

## How to interpret

- HRM separates by the **shape of the curve** (difference plot), not only by the peak Tm.
- The difference-plot peak (%) is an **uncalibrated, relative prediction** from a local
  melting model — use it to rank candidates, not as a calibrated value or a probability of
  success. Rule of thumb (generic heuristic): **>=25%** easy · **12–25%** likely ·
  **<12%** difficult. Confirm absolute behavior with uMelt and on your own instrument.
- Two samples with the **same genotype** (e.g. both heterozygous) **do not separate**.
- Markers with **STRONG** evidence do not depend on the "N" state (presence/absence);
  those that use "N" require bench confirmation of the allele.

## Limitations (honest)

- Indels are **inferred** from DArTseq presence/absence — **validate at the bench**
  (PCR + size) before finalizing the assay.
- The melting prediction is a **cooperative approximation** (no loop entropy à la
  Poland–Scheraga). For fine checking, use uMelt with the `amplicons.fasta`.
- Without a reference genome, primer yield is low (the tags are ~69 bp).

## License

MIT — see `LICENSE`. To cite, see `CITATION.cff`.
