"""Análise de HRM ciente do GENÓTIPO de cada amostra.

No HRM o que se compara é a curva de cada AMOSTRA conforme o genótipo dela naquele
marcador: A = homoduplex longo, B = homoduplex curto, AB = heterozigoto (homo + HETEROduplex,
derrete antes). Duas amostras com o mesmo genótipo têm a mesma curva e NÃO se separam.

Este módulo calcula, por marcador, a curva de cada estado (A/B/AB), a separação par a par
entre amostras (pico no difference plot) e constrói a árvore de decisão priorizando a
separação real de HRM. É genérico: funciona para qualquer conjunto de amostras.
"""
import itertools
from difflib import SequenceMatcher

import numpy as np

from .melt_model import helicity, het_curve

STRONG = {"A", "B", "AB"}


def indel_in_amplicon(amp_long, amp_short):
    """Posição e tamanho do indel dentro do amplicon (bases a mais no alelo longo)."""
    sm = SequenceMatcher(None, amp_long, amp_short, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "delete" and i2 > i1:
            return i1, i2 - i1
    return len(amp_long) // 2, 0


def _norm(theta):
    th = np.asarray(theta, dtype=float)
    lo, hi = th.min(), th.max()
    return (th - lo) / (hi - lo) if hi > lo else th * 0.0


def grid_of(cfg):
    return np.arange(cfg.melt_tmin, cfg.melt_tmax, cfg.melt_step)


def state_curves(amp_long, amp_short, cfg, grid):
    """Curvas normalizadas por estado {A,B,AB,N} na grade comum."""
    ipos, ilen = indel_in_amplicon(amp_long, amp_short)
    kw = dict(mono_mM=cfg.melt_mono, mg_mM=cfg.melt_mg, dntp_mM=cfg.melt_dntp,
              sigma=cfg.melt_sigma, tmin=cfg.melt_tmin, tmax=cfg.melt_tmax, step=cfg.melt_step)
    tA, thA, _, _ = helicity(amp_long, **kw)
    tB, thB, _, _ = helicity(amp_short, **kw)
    if ilen > 0:
        tH, thH, _, _ = het_curve(amp_long, ipos, ilen, mono_mM=cfg.melt_mono, mg_mM=cfg.melt_mg,
                                  dntp_mM=cfg.melt_dntp, sigma=cfg.melt_sigma,
                                  tmin=cfg.melt_tmin, tmax=cfg.melt_tmax, step=cfg.melt_step)
    else:
        tH, thH = tA, thA

    def itp(t, th):
        return _norm(np.interp(grid, t, th))

    cur = {"A": itp(tA, thA), "B": itp(tB, thB), "AB": itp(tH, thH)}
    cur["N"] = cur["B"]          # suposição: ausência de tag ~ alelo curto (sinalizado como não-forte)
    return cur, ipos, ilen


def diff_metrics(c1, c2):
    d = (np.asarray(c1) - np.asarray(c2)) * 100.0
    return round(float(np.max(np.abs(d))), 1), round(float(np.sqrt(np.mean(d ** 2))), 1)


def attach_curves(usable, cfg):
    """Adiciona 'curves' (por estado) a cada marcador usável (que tem primers)."""
    grid = grid_of(cfg)
    for r in usable:
        p = r["primers"]["pairs"][0]
        r["curves"], r["ipos_amp"], r["ilen_amp"] = state_curves(
            p["amplicon_long"], p["amplicon_short"], cfg, grid)
    return grid


def pair_separation(marker, sample_x, sample_y, labels_key="labels"):
    """Pico/RMS do difference plot entre duas amostras, pelo genótipo de cada uma."""
    sx = marker[labels_key].get(sample_x, "miss")
    sy = marker[labels_key].get(sample_y, "miss")
    if sx == "miss" or sy == "miss" or sx == sy:
        return 0.0, 0.0, True
    mx, rms = diff_metrics(marker["curves"][sx], marker["curves"][sy])
    return mx, rms, ("N" not in (sx, sy))


def best_marker_per_pair(usable, samples):
    """Para cada par de amostras, ranqueia os marcadores pela separação (FORTE primeiro)."""
    out = {}
    for x, y in itertools.combinations(samples, 2):
        rows = []
        for r in usable:
            mx, rms, conf = pair_separation(r, x, y)
            if mx > 0:
                rows.append(dict(marker=r, peak=mx, rms=rms, forte=conf,
                                 sx=r["labels"][x], sy=r["labels"][y]))
        rows.sort(key=lambda d: (not d["forte"], -d["peak"]))
        out[(x, y)] = rows
    return out


# ---------------- árvore de decisão ciente do genótipo ----------------

def _mkey(r):
    return (r["chrom"], r["pos"])


def _groups(marker, members):
    g = {}
    for a in members:
        s = marker["labels"].get(a, "miss")
        if s == "miss":
            return None
        g.setdefault(s, []).append(a)
    return g


def _node_sep(marker, members):
    """Menor separação entre membros de estados diferentes, e se é tudo FORTE."""
    wmin, forte = 999.0, True
    for a, b in itertools.combinations(members, 2):
        if marker["labels"][a] != marker["labels"][b]:
            sp, _, ft = pair_separation(marker, a, b)
            wmin = min(wmin, sp)
            forte = forte and ft
    return wmin, forte


def _best_subtree(members, used, pool):
    members = tuple(sorted(members))
    if len(members) <= 1:
        return dict(leaf=members[0] if members else None, markers=frozenset(),
                    depth=0, wmin=999.0, forte=True)
    best, bkey = None, None
    for r in pool:
        k = _mkey(r)
        if k in used:
            continue
        g = _groups(r, members)
        if not g or len(g) < 2:
            continue
        wmin, forte_here = _node_sep(r, members)
        subused = set(used) | {k}
        subs, ok = {}, True
        for s, mem in g.items():
            st = _best_subtree(mem, subused, pool)
            if st is None:
                ok = False
                break
            subs[s] = (mem, st)
            subused |= set(st["markers"])
        if not ok:
            continue
        mk = frozenset({k}) | frozenset().union(*[st["markers"] for _, st in subs.values()])
        depth = 1 + max(st["depth"] for _, st in subs.values())
        wtot = min([wmin] + [st["wmin"] for _, st in subs.values() if st["wmin"] < 999])
        forte = forte_here and all(st["forte"] for _, st in subs.values())
        cand = dict(node=r, groups=subs, markers=mk, depth=depth,
                    wmin=round(wtot, 1), forte=forte, nodesep=round(wmin, 1), leaf=None)
        key = (-wtot, not forte, depth)
        if bkey is None or key < bkey:
            bkey, best = key, cand
    return best


def build_trees(usable, samples, cfg, n_trees=3):
    """Constrói árvores (uma por raiz possível), ranqueia e escolhe até n_trees disjuntas."""
    pool = usable
    cands = []
    for r in pool:
        g = _groups(r, samples)
        if not g or len(g) < 2 or max(len(v) for v in g.values()) > 3:
            continue
        used = {_mkey(r)}
        subs, ok = {}, True
        for s, mem in g.items():
            st = _best_subtree(mem, used, pool)
            if st is None:
                ok = False
                break
            subs[s] = (mem, st)
            used |= set(st["markers"])
        if not ok:
            continue
        wmin, forte_here = _node_sep(r, samples)
        mk = frozenset({_mkey(r)}) | frozenset().union(*[st["markers"] for _, st in subs.values()])
        depth = 1 + max(st["depth"] for _, st in subs.values())
        wtot = min([wmin] + [st["wmin"] for _, st in subs.values() if st["wmin"] < 999])
        forte = forte_here and all(st["forte"] for _, st in subs.values())
        cands.append(dict(node=r, groups=subs, markers=mk, depth=depth,
                          wmin=round(wtot, 1), forte=forte, nodesep=round(wmin, 1)))
    uniq = {}
    for t in cands:
        uniq.setdefault(t["markers"], t)
    cands = sorted(uniq.values(), key=lambda t: (not t["forte"], -t["wmin"], t["depth"]))
    if not cands:
        return []
    chosen = [cands[0]]
    usedsets = [cands[0]["markers"]]
    roots = {_mkey(cands[0]["node"])}
    while len(chosen) < n_trees and len(chosen) < len(cands):
        best, bkey = None, None
        for t in cands:
            if t in chosen:
                continue
            ov = sum(len(t["markers"] & u) for u in usedsets)
            key = (_mkey(t["node"]) in roots, not t["forte"], ov, -t["wmin"], t["depth"])
            if bkey is None or key < bkey:
                bkey, best = key, t
        chosen.append(best)
        usedsets.append(best["markers"])
        roots.add(_mkey(best["node"]))
    return chosen
