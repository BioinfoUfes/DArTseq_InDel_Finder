"""Marker selection: ranking and minimal panel (greedy set cover).

Works on the characterized list. 'targets' are the samples of interest (all, or a group).
The "clean" separation uses only A/B/AB states. The genotype-aware decision tree lives in
genotype.py.
"""
import itertools


def _clean_map(o, targets):
    return {c: o['labels'][c] for c in targets if o['labels'][c] in ('A', 'B', 'AB')}


def _pairset(o, targets, pair_index):
    cm = _clean_map(o, targets)
    s = set()
    for idx, (c1, c2) in pair_index.items():
        if c1 in cm and c2 in cm and cm[c1] != cm[c2]:
            s.add(idx)
    return s


def rank(out, cfg):
    def key(o):
        sil = 0 if o['src'] == 'SIL-SIL' else (1 if 'SIL' in o['src'] else 2)
        good_size = 0 if 3 <= o['indel'] <= 15 else 1
        return (sil, -o['sep_clean'], o['hp'], good_size, o['incons_n'], -o['ident'])
    return sorted(out, key=key)


def minimal_panel(out, targets, cfg):
    """Smallest set of InDels (greedy) that separates every pair of 'targets'."""
    pairs = list(itertools.combinations(targets, 2))
    pair_index = {i: p for i, p in enumerate(pairs)}
    pool = [o for o in out
            if (not cfg.only_sil_sil or o['src'] == 'SIL-SIL')
            and (o['hp'] is False or not cfg.exclude_homopolymer)
            and o['indel'] >= cfg.indel_min
            and o['incons_n'] <= max(2, len(targets) // 3 + 1)]
    psets = [(_pairset(o, targets, pair_index), o) for o in pool]
    cover_any = set().union(*[p for p, _ in psets]) if psets else set()
    all_idx = set(pair_index)
    covered = set()
    panel = []
    used = [False] * len(psets)
    while True:
        best = -1; bi = -1
        for k, (ps, _) in enumerate(psets):
            if used[k]:
                continue
            g = len(ps - covered)
            if g > best:
                best = g; bi = k
        if best <= 0:
            break
        used[bi] = True
        panel.append((psets[bi][1], best))
        covered |= psets[bi][0]
    # unresolved = every pair the final panel does not separate, INCLUDING pairs that no
    # candidate marker can separate (resolvable < n_pairs reveals those).
    unresolved = sorted(all_idx - covered)
    return dict(panel=panel, n_pairs=len(pairs),
                resolved=len(covered),
                resolvable=len(cover_any),
                unresolved=[pair_index[i] for i in unresolved])
