"""Orchestrator: from DArTseq reports to HRM decision trees with primers.

Generic: works for any species — just provide the DArTseq reports and one (or more)
reference genome FASTA(s). Markers are located in the genome by the tag sequence itself
(remap), so it works for any assembly.
"""
import itertools

from .io_report import load_reports
from .indel_finder import find_indels
from .characterize import characterize
from .select import rank, minimal_panel
from .region import build_regions, _load_fasta_index
from .primers import design_for_region
from .insilico import screen_markers
from .genotype import attach_curves, build_trees, best_marker_per_pair


def _usable_filter(region):
    """Keep markers with a specific primer (single copy in the genome)."""
    if not (region["primers"]["ok"] and region["primers"]["pairs"]):
        return False
    if not region["region_source"].startswith("genome"):
        return False
    gh = region.get("genome_hits")
    if isinstance(gh, dict) and gh and max(gh.values()) > 1:
        return False
    return True


def _candidate_pool(out, samples, cfg, cap=80):
    """Markers that separate at least one pair of 'samples' (or can serve as a root)."""
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
    # largest indel first (proxy for larger separation), then highest discrimination
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


def _design_and_keep(markers, cfg, genome_index, log=lambda *_: None):
    regions = build_regions(markers, cfg, genome_index)
    usable = []
    for r in regions:
        r["primers"] = design_for_region(r, cfg)
        if _usable_filter(r):
            usable.append(r)
    if usable and genome_index is not None:
        log(f"  in-silico PCR screen on {len(usable)} markers...")
        screen_markers(usable, genome_index, cfg)
        if cfg.insilico_filter:
            kept = [r for r in usable if r.get("insilico", {}).get("specific")]
            log(f"  in-silico specific: {len(kept)}/{len(usable)} (filter ON)")
            usable = kept
    return usable


def run(cfg, samples=None, mode="tree", n_trees=3, log=print):
    """Run the pipeline. 'samples' = list of sample codes (None = all).
    mode='tree' (decision tree) or 'panel' (minimal panel). Returns a results dict."""
    log("Reading DArTseq reports...")
    data = load_reports(cfg.reports_dir)
    all_samples = data["cult"]
    samples = samples or all_samples
    log(f"  samples in report: {len(all_samples)} | target: {len(samples)}")

    log("Detecting indel candidates...")
    cand = find_indels(data["snp"] + data["sil"], cfg)
    out, _ = characterize(cand, data["snp"] + data["sil"], all_samples)
    out = rank(out, cfg)
    log(f"  indel candidates: {len(out)}")

    genome_index = None
    if cfg.genome_fna:
        log("Loading reference genome(s)...")
        genome_index = _load_fasta_index(cfg.genome_fna)
        log(f"  sequences in index: {len(genome_index)}")
    else:
        log("  WARNING: no genome — using tags only (~69 bp); low primer yield.")

    result = dict(samples=samples, all_samples=all_samples, mode=mode, out=out)

    if mode == "panel":
        mp = minimal_panel(out, samples, cfg)
        markers = [o for o, _ in mp["panel"]]
        log(f"Designing primers for {len(markers)} panel markers...")
        usable = _design_and_keep(markers, cfg, genome_index, log) if genome_index else []
        if usable:
            attach_curves(usable, cfg)
        result.update(panel=mp, usable=usable)
        return result

    # tree mode
    pool = _candidate_pool(out, samples, cfg)
    log(f"Designing primers (genome) for {len(pool)} candidates...")
    usable = _design_and_keep(pool, cfg, genome_index, log) if genome_index else []
    log(f"  usable markers (specific primer): {len(usable)}")
    if not usable:
        log("  No usable markers — check the reference genome.")
        result.update(usable=[], trees=[], per_pair={})
        return result
    attach_curves(usable, cfg)
    log("Computing genotype-based HRM separation and building trees...")
    trees = build_trees(usable, samples, cfg, n_trees=n_trees)
    per_pair = best_marker_per_pair(usable, samples)
    result.update(usable=usable, trees=trees, per_pair=per_pair)
    return result
