"""Central parameters of the InDel-to-HRM pipeline.

Edit the quality thresholds, amplicon size and primer criteria here. No path is hard-coded
anywhere outside this file.
"""
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Config:
    # ---- input ----
    reports_dir: Path = Path("reports")          # Report_*_SNP_mapping_*.csv and *_SilicoDArT_*.csv
    genome_fna: Path | None = None               # reference genome FASTA(s); None = tags only
    out_dir: Path = Path("indel_hrm_output")

    # ---- marker quality (DArTseq) ----
    min_rep: float = 0.95        # minimum RepAvg / Reproducibility
    min_callrate: float = 0.80   # minimum CallRate
    min_pic: float = 0.0         # minimum PIC (0 = no filter here)

    # ---- InDel detection ----
    coloc_window: int = 60       # tags <= N bp apart in the genome are candidates for the same locus
    indel_min: int = 3           # minimum InDel size (bp)
    indel_max: int = 45          # maximum size (larger InDels = more distinguishable by HRM)
    min_ident: float = 0.82      # minimum identity in the aligned flanks
    min_anchor: int = 20         # minimum number of matched bases (anchor)
    min_exact_anchor: int = 0    # minimum EXACT (identical) run flanking the gap, each side;
                                 # 0 = off (default). Set 15-25 for stricter, repetitive genomes
    max_gap_runs: int = 2        # maximum number of gap runs (parsimony)
    exclude_homopolymer: bool = True

    # ---- selection ----
    require_triplicate_consistency: bool = True  # require agreeing replicates
    only_sil_sil: bool = True    # use only SilicoDArT x SilicoDArT pairs (best InDel signal)
    insilico_filter: bool = False  # if True, drop markers whose primer pair is not
                                   # in-silico specific (>1 product in a single haplotype);
                                   # default False = report the count, do not filter

    # ---- flanks for the primer template ----
    flank: int = 180             # bp on each side of the InDel taken from the genome (when provided)
    remap: bool = True           # locate the tag in the genome by sequence (ignores ChromPos);
                                 # required when using a genome other than the one DArT aligned to
    remap_seed: int = 25         # anchor size (bp next to the InDel) used in the search
    remap_max_hits: int = 1      # max number of genome hits to accept (specificity);
                                 # tags matching more places are flagged as non-specific

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
    n_primer_pairs: int = 3      # primer pairs returned per marker

    # ---- primer Tm (Primer3) ----
    mono_conc: float = 50.0      # [Na+] mM for the primer Tm calculation
    div_conc: float = 3.0        # [Mg2+] mM
    dntp_conc: float = 0.8       # mM
    dna_conc: float = 250.0      # nM

    # ---- HRM melting-curve prediction (Zimm-Bragg/Ising model) ----
    # Set these to YOUR HRM kit / master-mix conditions.
    melt_mono: float = 20.0      # [Mono+] mM of the HRM master mix
    melt_mg: float = 3.0         # [Mg2+] mM of the HRM master mix
    melt_dntp: float = 0.8       # [dNTP] mM
    melt_sigma: float = 1e-4     # cooperativity (Zimm-Bragg)
    melt_tmin: float = 45.0      # temperature range of the simulation (C)
    melt_tmax: float = 98.0
    melt_step: float = 0.2       # resolution (C)

    # ---- optional input ----
    metadata_csv: Path | None = None  # CSV with columns 'code,name' to label the samples
    # HRM interpretation thresholds (difference-plot peak, %)
    hrm_easy: float = 25.0       # >= => easy separation
    hrm_likely: float = 12.0     # >= => likely; below => difficult
