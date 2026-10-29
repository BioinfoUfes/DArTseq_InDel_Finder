# Changelog

## v2.0.0 — 2026-10-29

Methodological and software-quality upgrade. The command-line interface and the overall
workflow are unchanged; the decision-tree result on the reference dataset is preserved.

### Added
- **Affine-gap DNA aligner** (Biopython `PairwiseAligner`, free terminal gaps) as the
  aligner of record for InDel calling, with explicit match/mismatch/gap scoring. A cheap
  pre-gate keeps it fast and a difflib fallback keeps the module working without Biopython.
  An exact-flanking-anchor metric is recorded per candidate (optional stricter filter via
  `--min-exact-anchor`).
- **Exact-match in-silico PCR screen** against the reference, counting predicted products
  per haplotype (homeolog amplifying once per haplotype is expected; >1 in one haplotype
  flags a non-specific pair). Reported in the outputs; optional hard filter via
  `--insilico-filter`.
- **Machine-readable outputs** alongside the XLSX: `candidates.tsv`, `primers.tsv`,
  `per_pair.tsv`, and `run_summary.json` (parameters, counts and chosen tree).
- **Unit tests** (`tests/`), a **GitHub Actions** workflow, and **`pyproject.toml`**
  packaging. Installable console entry points: `dartseq-indel-finder` and `python -m indel_hrm`.

### Changed
- **Honest handling of `N` (double null).** `N` is no longer aliased to allele `B` when
  predicting curves; it has no predicted curve, and any pair that relies on `N` is reported
  as "uses N (confirm)", never as STRONG.
- **Honest labeling of HRM predictions.** Difference-plot peaks are labeled as uncalibrated,
  relative predictions from a local melting model — not instrument-calibrated values and not
  a probability of assay success. The easy/likely/difficult thresholds are documented as
  generic heuristics; uMelt and bench validation remain the references for absolute behavior.
- Repository fully in English (code comments, interactive prompts, reports and outputs).

### Preserved
- Haplotype-aware tag specificity (homeolog vs paralog) — better suited to phased assemblies
  than a blunt unique-mapping gate.
- Genotype-aware decision trees and the greedy minimal panel, with unresolved pairs reported.
- Same authors, name and citation lineage as `DArTseq_SNP_Finder`.

## v1.0.0 — 2026-10-29
- Initial deposit (Zenodo concept DOI 10.5281/zenodo.23775496).
