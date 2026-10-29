# Exemplos

`amostras_exemplo.csv` — modelo de CSV de nomes (colunas `code,name`). Os códigos
devem bater com os nomes de amostra do relatório DArTseq (ex.: `genotype 7011` → `7011`).

## Exemplo de execução (grupo de 4 acessos de Ponkan)

```bash
python ../run.py \
  --reports ../reports \
  --genome ../genomes/hapA.fa.gz ../genomes/hapB.fa.gz \
  --metadata amostras_exemplo.csv \
  --samples 7011,7041,7115,7167 \
  --mode tree --melt-mg 3 --melt-mono 20 \
  --out saida_ponkan
```

Coloque os relatórios DArTseq em `reports/` e os genomas em `genomes/`
(ambos são ignorados pelo git — veja `.gitignore`).
