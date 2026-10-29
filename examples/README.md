# Examples

`example_samples.csv` — template name CSV (columns `code,name`). The codes must match
the sample names in the DArTseq report (e.g. `genotype 7011` → `7011`).

## Example run (group of 4 Ponkan accessions)

```bash
python ../run.py \
  --reports ../reports \
  --genome ../genomes/hapA.fa.gz ../genomes/hapB.fa.gz \
  --metadata example_samples.csv \
  --samples 7011,7041,7115,7167 \
  --mode tree --melt-mg 3 --melt-mono 20 \
  --out output_ponkan
```

Put the DArTseq reports in `reports/` and the genomes in `genomes/`
(both are ignored by git — see `.gitignore`).
