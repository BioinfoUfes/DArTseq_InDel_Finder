"""Parâmetros centrais da pipeline de indels para HRM.

Edite aqui os limiares de qualidade, tamanho de amplicon e critérios de primer.
Nenhum caminho é fixo em código fora deste arquivo.
"""
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Config:
    # ---- entrada ----
    reports_dir: Path = Path("reports")          # Report_*_SNP_mapping_*.csv e *_SilicoDArT_*.csv
    genome_fna: Path | None = None               # GCF_000317415.1 (Csi_valencia_1.0); None = usa só as tags
    out_dir: Path = Path("indel_hrm_output")

    # ---- qualidade do marcador (DArTseq) ----
    min_rep: float = 0.95        # RepAvg / Reproducibility mínimo
    min_callrate: float = 0.80   # CallRate mínimo
    min_pic: float = 0.0         # PIC mínimo (0 = não filtra aqui)

    # ---- detecção de indel ----
    coloc_window: int = 60       # tags a <= N pb no genoma são candidatas ao mesmo locus
    indel_min: int = 3           # tamanho mínimo do indel (pb)
    indel_max: int = 45          # tamanho máximo (indels maiores = mais distinguíveis no HRM)
    min_ident: float = 0.82      # identidade mínima nos flancos alinhados
    min_anchor: int = 20         # nº mínimo de bases pareadas (âncora)
    max_gap_runs: int = 2        # nº máximo de corridas de gap (parcimônia)
    exclude_homopolymer: bool = True

    # ---- seleção ----
    require_triplicate_consistency: bool = True  # exige réplicas concordantes
    only_sil_sil: bool = True    # usar só pares SilicoDArT×SilicoDArT (melhor sinal de indel)

    # ---- flancos para o template de primer ----
    flank: int = 180             # pb de cada lado do indel extraídos do genoma (quando houver .fna)
    remap: bool = True           # localizar a tag no genoma pela sequência (ignora ChromPos);
                                 # obrigatório ao usar um genoma diferente do Csi_valencia_1.0
    remap_seed: int = 25         # tamanho da âncora (pb a montante do indel) usada na busca
    remap_max_hits: int = 1      # nº máx. de locais no genoma p/ aceitar (especificidade);
                                 # tags que batem em mais lugares são marcadas como não-específicas

    # ---- Primer3 / HRM ----
    hrm_product_min: int = 50
    hrm_product_max: int = 120
    hrm_product_opt: int = 90
    primer_opt_size: int = 20
    primer_min_size: int = 17
    primer_max_size: int = 25
    primer_opt_tm: float = 60.0
    primer_min_tm: float = 57.0
    primer_max_tm: float = 63.0
    primer_min_gc: float = 30.0
    primer_max_gc: float = 70.0
    primer_max_pair_tm_diff: float = 1.5
    n_primer_pairs: int = 3      # pares retornados por marcador

    # ---- Tm dos primers (Primer3) ----
    mono_conc: float = 50.0      # [Na+] mM para o cálculo de Tm do primer
    div_conc: float = 3.0        # [Mg2+] mM
    dntp_conc: float = 0.8       # mM
    dna_conc: float = 250.0      # nM

    # ---- predição da curva de HRM (modelo Zimm-Bragg/Ising) ----
    # Use aqui as condições do SEU kit/master mix de HRM.
    melt_mono: float = 20.0      # [Mono+] mM do master mix de HRM
    melt_mg: float = 3.0         # [Mg2+] mM do master mix de HRM
    melt_dntp: float = 0.8       # [dNTP] mM
    melt_sigma: float = 1e-4     # cooperatividade (Zimm-Bragg)
    melt_tmin: float = 45.0      # faixa de temperatura da simulação (°C)
    melt_tmax: float = 98.0
    melt_step: float = 0.2       # resolução (°C)

    # ---- entrada opcional ----
    metadata_csv: Path | None = None  # CSV com colunas 'code,name' p/ rotular as amostras
    # limiares de interpretação do HRM (pico no difference plot, %)
    hrm_easy: float = 25.0       # >= => separação fácil
    hrm_likely: float = 12.0     # >= => provável; abaixo => difícil
