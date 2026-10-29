"""Perguntas interativas para montar a configuração da pipeline."""
from pathlib import Path

from .config import Config
from .io_report import load_reports


def _ask(prompt, default=None, cast=str):
    sfx = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{prompt}{sfx}: ").strip()
        if not raw and default is not None:
            return default
        if not raw:
            print("  (resposta obrigatória)")
            continue
        try:
            return cast(raw)
        except Exception:
            print("  (valor inválido, tente de novo)")


def _ask_path(prompt, default=None, must_exist=True):
    while True:
        p = _ask(prompt, default)
        path = Path(p).expanduser()
        if not must_exist or path.exists():
            return path
        print(f"  (caminho não encontrado: {path})")


def interactive_config():
    print("=" * 64)
    print(" indel-hrm-finder — configuração interativa")
    print("=" * 64)
    reports = _ask_path("Pasta com os relatórios DArTseq (SNP_mapping + SilicoDArT)",
                        "reports")
    # genomas
    print("\nGenoma(s) de referência da sua espécie (FASTA .fa/.fa.gz).")
    print("Pode informar vários separados por vírgula (ex.: dois haplótipos). Enter = nenhum.")
    graw = input("Arquivo(s) de genoma: ").strip()
    genome = [Path(g.strip()).expanduser() for g in graw.split(",") if g.strip()] or None
    if genome:
        for g in genome:
            if not g.exists():
                print(f"  AVISO: não encontrei {g}")

    meta = input("\nCSV de nomes das amostras (colunas code,name) — Enter para pular: ").strip()
    metadata = Path(meta).expanduser() if meta else None

    mode = _ask("\nModo: 'tree' (árvore p/ um grupo) ou 'panel' (painel p/ todas)", "tree")
    while mode not in ("tree", "panel"):
        mode = _ask("Digite 'tree' ou 'panel'", "tree")

    # amostras disponíveis
    try:
        data = load_reports(reports)
        codes = data["cult"]
        print(f"\nAmostras no relatório ({len(codes)}): {', '.join(codes)}")
    except Exception as e:
        print(f"  (não consegui listar amostras agora: {e})")
        codes = []

    samples = None
    if mode == "tree":
        sraw = input("Códigos das amostras do grupo (separados por vírgula) ou 'all': ").strip()
        if sraw and sraw.lower() != "all":
            samples = [s.strip() for s in sraw.split(",") if s.strip()]

    print("\nCondições do seu kit/master mix de HRM (Enter = padrão):")
    mg = _ask("  [Mg2+] mM", 3.0, float)
    mono = _ask("  [Mono+] mM", 20.0, float)

    indel_max = _ask("\nTamanho máximo de indel a buscar (pb)", 45, int)
    flank = _ask("Flanco extraído do genoma p/ desenhar primer (pb de cada lado)", 180, int)
    out = _ask("\nPasta de saída", "indel_hrm_output")

    cfg = Config(reports_dir=reports, genome_fna=genome, out_dir=Path(out),
                 metadata_csv=metadata, indel_max=indel_max, flank=flank,
                 melt_mg=mg, melt_mono=mono)
    return cfg, samples, mode
