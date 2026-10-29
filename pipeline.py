"""Orquestrador: dos relatórios DArTseq às árvores de decisão de HRM com primers.

Genérico: funciona para qualquer espécie, bastando fornecer os relatórios DArTseq e um
(ou mais) genoma(s) de referência FASTA. A localização dos marcadores no genoma é feita
pela própria sequência da tag (remap), então serve para qualquer montagem.
"""
import itertools

from .io_report import load_reports
from .indel_finder import find_indels
from .characterize import characterize
from .select import rank, minimal_panel
from .region import build_regions, _load_fasta_index
from .primers import design_for_region
from .genotype import attach_curves, build_trees, best_marker_per_pair


def _usable_filter(region):
    """Mantém marcadores com primer específico (cópia única no genoma)."""
    if not (region["primers"]["ok"] and region["primers"]["pairs"]):
        return False
    if not region["region_source"].startswith("genome"):
        return False
    gh = region.get("genome_hits")
    if isinstance(gh, dict) and gh and max(gh.values()) > 1:
        return False
    return True


def _candidate_pool(out, samples, cfg, cap=80):
    """Marcadores que separam ao menos um par dos 'samples' (ou servem de raiz)."""
    pairs = list(itertools.combinations(samples, 2))

    def st_of(o):
        return {a: o["labels"][a] for a in samples
                if o["labels"][a] != "miss" and not o["incons"][a]}

    pool = []
    for o in out:
        if o["src"] != "SIL-SIL" or o["hp"] or o["indel"] < cfg.indel_min:
            continue
        st = st_of(o)
        seps = any(a in st and b in st and st[a] != st[b] for a, b in pairs)
        valid_root = (len(st) == len(samples)
                      and len({st[a] for a in samples}) >= 2)
        if seps or valid_root:
            pool.append(o)
    # indel grande primeiro (proxy de maior separação), depois maior discriminação
    pool.sort(key=lambda o: (-o["indel"], -o["sep_clean"]))
    seen, uniq = set(), []
    for o in pool:
        k = (o["chrom"], o["pos"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(o)
        if len(uniq) >= cap:
            break
    return uniq


def _design_and_keep(markers, cfg, genome_index):
    regions = build_regions(markers, cfg, genome_index)
    usable = []
    for r in regions:
        r["primers"] = design_for_region(r, cfg)
        if _usable_filter(r):
            usable.append(r)
    return usable


def run(cfg, samples=None, mode="tree", n_trees=3, log=print):
    """Executa a pipeline. 'samples' = lista de códigos de amostra (None = todas).
    mode='tree' (árvore de decisão) ou 'panel' (painel mínimo). Retorna dict de resultados."""
    log("Lendo relatórios DArTseq...")
    data = load_reports(cfg.reports_dir)
    all_samples = data["cult"]
    samples = samples or all_samples
    log(f"  amostras no relatório: {len(all_samples)} | alvo: {len(samples)}")

    log("Detectando candidatos a indel...")
    cand = find_indels(data["snp"] + data["sil"], cfg)
    out, _ = characterize(cand, data["snp"] + data["sil"], all_samples)
    out = rank(out, cfg)
    log(f"  candidatos a indel: {len(out)}")

    genome_index = None
    if cfg.genome_fna:
        log("Carregando genoma(s) de referência...")
        genome_index = _load_fasta_index(cfg.genome_fna)
        log(f"  sequências no índice: {len(genome_index)}")
    else:
        log("  AVISO: sem genoma — usando só as tags (~69 pb); rendimento de primer baixo.")

    result = dict(samples=samples, all_samples=all_samples, mode=mode, out=out)

    if mode == "panel":
        mp = minimal_panel(out, samples, cfg)
        markers = [o for o, _ in mp["panel"]]
        log(f"Desenhando primers para {len(markers)} marcadores do painel...")
        usable = _design_and_keep(markers, cfg, genome_index) if genome_index else []
        if usable:
            attach_curves(usable, cfg)
        result.update(panel=mp, usable=usable)
        return result

    # modo árvore
    pool = _candidate_pool(out, samples, cfg)
    log(f"Desenhando primers (genoma) para {len(pool)} candidatos...")
    usable = _design_and_keep(pool, cfg, genome_index) if genome_index else []
    log(f"  marcadores usáveis (primer específico): {len(usable)}")
    if not usable:
        log("  Nenhum marcador usável — verifique o genoma de referência.")
        result.update(usable=[], trees=[], per_pair={})
        return result
    attach_curves(usable, cfg)
    log("Calculando separação de HRM por genótipo e montando árvores...")
    trees = build_trees(usable, samples, cfg, n_trees=n_trees)
    per_pair = best_marker_per_pair(usable, samples)
    result.update(usable=usable, trees=trees, per_pair=per_pair)
    return result
