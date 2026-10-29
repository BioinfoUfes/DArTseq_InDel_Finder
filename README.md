# **DArTseq_InDel_Finder**
---
### InDels that Differentiate Samples (Cultivars / Species)

A pipeline to detect, filter, and select **InDels** obtained from DArTseq that
differentiate groups of samples (e.g., cultivars or species). Beyond detecting the
InDels, it designs HRM primers and predicts, per genotype, the High Resolution Melting
separation — so the output is a set of prioritized, ready-to-test HRM assays.

> É a evolução (InDel + HRM) do [`DArTseq_SNP_Finder`](https://github.com/BioinfoUfes/DArTseq_SNP_Finder),
> desenvolvida pelo mesmo grupo (BioInfo/UFES).

A partir dos relatórios DArTseq e de um (ou mais) genoma(s) de referência da sua
espécie, a ferramenta:

1. **infere indels** comparando tags que caem no mesmo locus do genoma (o DArTseq-SNP
   só reporta substituições; os indels são inferidos);
2. **mapeia e checa especificidade** de cada tag no genoma que você fornecer (funciona
   com qualquer montagem, inclusive haplótipos separados);
3. **desenha primers** de HRM (Primer3) com amplicon curto, flanqueando o indel;
4. **prevê a curva de melting de cada amostra pelo genótipo** (A = homoduplex longo,
   B = homoduplex curto, AB = heterozigoto com heteroduplex), via um modelo local
   (Zimm-Bragg/Ising com parâmetros de vizinho-mais-próximo) — sem depender de sites;
5. **monta árvores de decisão** que identificam as amostras no menor número de análises,
   priorizando a separação real de HRM (pico no *difference plot*).

É **genérico**: serve para qualquer espécie, bastando trocar os relatórios e o genoma.

---

## Instalação

Com conda (recomendado):

```bash
conda env create -f environment.yml
conda activate indel-hrm-finder
```

Ou com pip:

```bash
pip install -r requirements.txt
```

Dependências: `primer3-py`, `biopython`, `numpy`, `openpyxl`, `matplotlib`.
Não é necessário BLAST nem instalar nada além disso.

## Entradas

- **Relatórios DArTseq** numa pasta (ex.: `reports/`): um `*SNP_mapping*.csv` **e** um
  `*SilicoDArT*.csv` do seu pedido DArT.
- **Genoma(s) de referência** da sua espécie, em FASTA (`.fa`/`.fa.gz`). Pode informar
  mais de um arquivo (ex.: dois haplótipos de uma montagem faseada). A localização é
  feita pela própria sequência da tag, então **qualquer montagem serve**.
- *(opcional)* **CSV de nomes** das amostras, com colunas `code,name`, para que os
  relatórios saiam com nomes legíveis em vez dos códigos (veja `examples/`).

## Uso

### Modo interativo (faz as perguntas)

```bash
python run.py
```

Ele pergunta: a pasta dos relatórios, o(s) genoma(s), o CSV de nomes (opcional), o modo
(árvore ou painel), quais amostras comparar, as condições de HRM do seu kit
(Mg²⁺ e Mono⁺) e a pasta de saída.

### Modo linha de comando

```bash
python run.py \
  --reports reports \
  --genome genomaHapA.fa.gz genomaHapB.fa.gz \
  --metadata examples/amostras_exemplo.csv \
  --samples 7011,7041,7115,7167 \
  --mode tree \
  --melt-mg 3 --melt-mono 20 \
  --out saida
```

- `--mode tree` monta a **árvore de decisão** para o grupo de amostras indicado em
  `--samples` (omita para usar todas).
- `--mode panel` monta o **painel mínimo** que distingue todas as amostras.
- `--melt-mg` / `--melt-mono`: use as concentrações do **seu** master mix de HRM.

## Saídas (na pasta `--out`)

- `HRM_marcadores.xlsx`
  - **Separacao_por_par** — para cada par de amostras, os melhores marcadores e o pico
    de separação no *difference plot* (sua ferramenta para escolher primers).
  - **Arvores** — as árvores de decisão (melhor primeiro), com o fluxo passo a passo.
  - **Primers** — primers (F/R, Tm), amplicons e especificidade.
  - **Metodo** — parâmetros e notas.
- `curvas/` — uma imagem por marcador: a curva de **cada amostra pelo genótipo** +
  o *difference plot*.
- `amplicons.fasta` — alelos longo/curto de cada marcador, para conferência no uMelt.

## Como interpretar

- O HRM separa pela **forma da curva** (difference plot), não só pela Tm do pico.
- Regra prática do pico no difference plot: **≥25%** fácil · **12–25%** provável ·
  **<12%** difícil.
- Duas amostras com o **mesmo genótipo** (ex.: ambas heterozigotas) **não se separam**.
- Marcadores de evidência **FORTE** não dependem do estado "N" (presença/ausência);
  os que usam "N" exigem confirmação de bancada do alelo.

## Limitações (honestas)

- Os indels são **inferidos** de presença/ausência do DArTseq — **valide na bancada**
  (PCR + tamanho) antes de fechar o ensaio.
- A predição de melting é uma **aproximação cooperativa** (sem entropia de loop à la
  Poland–Scheraga). Para conferência fina, use o uMelt com o `amplicons.fasta`.
- Sem genoma de referência, o rendimento de primer é baixo (as tags têm ~69 pb).

## Licença

MIT — veja `LICENSE`. Para citar, veja `CITATION.cff`.
