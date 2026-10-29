"""Seleção de marcadores: ranking, painel mínimo (cobertura gulosa) e árvore de decisão.

Trabalha sobre a lista caracterizada. 'targets' são as cultivares de interesse
(todas, ou um grupo como os 4 Ponkans). A separação "limpa" usa só estados A/B/AB.
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
    """Menor conjunto de indels (guloso) que separa todos os pares de 'targets'."""
    pairs = list(itertools.combinations(targets, 2))
    pair_index = {i: p for i, p in enumerate(pairs)}
    pool = [o for o in out
            if (not cfg.only_sil_sil or o['src'] == 'SIL-SIL')
            and (o['hp'] is False or not cfg.exclude_homopolymer)
            and o['indel'] >= cfg.indel_min
            and o['incons_n'] <= max(2, len(targets) // 3 + 1)]
    psets = [(_pairset(o, targets, pair_index), o) for o in pool]
    cover_any = set().union(*[p for p, _ in psets]) if psets else set()
    remaining = set(cover_any)
    panel = []
    used = [False] * len(psets)
    while remaining:
        best = -1; bi = -1
        for k, (ps, _) in enumerate(psets):
            if used[k]:
                continue
            g = len(ps & remaining)
            if g > best:
                best = g; bi = k
        if best <= 0:
            break
        used[bi] = True
        panel.append((psets[bi][1], best))
        remaining -= psets[bi][0]
    return dict(panel=panel, n_pairs=len(pairs),
                resolved=len(cover_any) - len(remaining),
                resolvable=len(cover_any),
                unresolved=[pair_index[i] for i in remaining])


def decision_tree(out, targets, cfg):
    """Árvore que identifica qualquer 'target' no menor nº de análises.
    Raiz: marcador consistente nos alvos que deixa <=2 por galho. Nível 2: discriminador do par."""
    def consistent_all(o):
        return all((not o['incons'][c]) and o['labels'][c] != 'miss' for c in targets)

    def strong(o, members):
        return all(o['labels'][c] not in ('N', 'miss') for c in members)

    pool = [o for o in out if (not cfg.only_sil_sil or o['src'] == 'SIL-SIL')
            and not (cfg.exclude_homopolymer and o['hp']) and o['indel'] >= cfg.indel_min]

    roots = []
    for o in pool:
        if not consistent_all(o):
            continue
        groups = {}
        for c in targets:
            groups.setdefault(o['labels'][c], []).append(c)
        if max(len(v) for v in groups.values()) > 2 or len(groups) < 2:
            continue
        roots.append((o, groups))
    roots.sort(key=lambda og: (not strong(og[0], targets), -len(og[1]),
                               -min(og[0]['crA'], og[0]['crB']), -og[0]['indel']))

    def disc(c1, c2):
        opts = []
        for o in pool:
            if (not o['incons'][c1] and not o['incons'][c2]
                    and o['labels'][c1] != 'miss' and o['labels'][c2] != 'miss'
                    and o['labels'][c1] != o['labels'][c2]):
                st = o['labels'][c1] != 'N' and o['labels'][c2] != 'N'
                opts.append((o, st))
        opts.sort(key=lambda os_: (not os_[1], -min(os_[0]['crA'], os_[0]['crB']), -os_[0]['indel']))
        return opts

    best = None
    for o, groups in roots[:80]:
        subtests = []; ok = True
        for st, members in groups.items():
            if len(members) == 2:
                d = disc(members[0], members[1])
                if not d:
                    ok = False; break
                subtests.append((tuple(members), d[0]))
        if not ok:
            continue
        allstrong = strong(o, targets) and all(s for _, (_, s) in subtests)
        score = (allstrong, strong(o, targets), len(groups), min(o['crA'], o['crB']))
        if best is None or score > best['score']:
            best = dict(root=o, groups=groups, subtests=subtests, score=score)
    return best
