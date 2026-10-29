"""Detecção de candidatos a indel a partir de tags DArTseq co-localizadas.

Como o DArTseq-SNP só reporta substituições, os indels são inferidos: agrupamos
tags que o BLAST do DArT posicionou no mesmo locus (mesmo contig, posições a <=
coloc_window pb) e alinhamos cada par; um gap INTERNO (não terminal) de tamanho
indel_min..indel_max, com flancos de alta identidade, é um candidato a indel.
"""
from collections import defaultdict
from difflib import SequenceMatcher


def _indel_check(a, b):
    sm = SequenceMatcher(None, a, b, autojunk=False)
    ops = sm.get_opcodes()
    match = mism = internal = 0
    gaps = []
    n = len(ops)
    for oi, (tag, i1, i2, j1, j2) in enumerate(ops):
        if tag == 'equal':
            match += (i2 - i1)
        elif tag == 'replace':
            mism += max(i2 - i1, j2 - j1)
        elif tag in ('delete', 'insert'):
            g = (i2 - i1) + (j2 - j1)
            if oi not in (0, n - 1):          # ignora overhang terminal (trimagem)
                internal += g
                gaps.append((tag, g, i1, i2, j1, j2))
    aligned = match + mism
    ident = match / aligned if aligned else 0.0
    return match, mism, internal, gaps, ident


def find_indels(markers, cfg):
    by = defaultdict(list)
    for k, m in enumerate(markers):
        by[m['chrom']].append(k)
    cand = []
    for chrom, idxs in by.items():
        idxs.sort(key=lambda k: markers[k]['pos'])
        for ai in range(len(idxs)):
            ka = idxs[ai]; ma = markers[ka]
            if ma['rep'] < cfg.min_rep or ma['cr'] < cfg.min_callrate:
                continue
            pa = ma['pos']
            for bi in range(ai + 1, len(idxs)):
                kb = idxs[bi]; mb = markers[kb]
                if mb['pos'] - pa > cfg.coloc_window:
                    break
                if mb['rep'] < cfg.min_rep or mb['cr'] < cfg.min_callrate:
                    continue
                if cfg.only_sil_sil and not (ma['src'] == 'SIL' and mb['src'] == 'SIL'):
                    continue
                match, mism, internal, gaps, ident = _indel_check(ma['seq'], mb['seq'])
                if internal < cfg.indel_min or internal > cfg.indel_max:
                    continue
                if ident < cfg.min_ident or match < cfg.min_anchor:
                    continue
                if len(gaps) > cfg.max_gap_runs:
                    continue
                cand.append(dict(ka=ka, kb=kb, match=match, mism=mism,
                                 indel=internal, gaps=gaps, ident=round(ident, 3)))
    return cand


def gap_segment_and_homopolymer(seqA, seqB, gaps):
    """Retorna (sequência do indel, flag homopolímero/STR)."""
    hp = False
    segs = []
    for tag, g, i1, i2, j1, j2 in gaps:
        if tag == 'delete':
            seg = seqA[i1:i2]; pos = i1; ref = seqA
        else:
            seg = seqB[j1:j2]; pos = j1; ref = seqB
        segs.append(seg)
        if seg and len(set(seg)) == 1:                  # homopolímero
            base = seg[0]; run = len(seg); k = pos - 1
            while k >= 0 and ref[k] == base:
                run += 1; k -= 1
            k = pos + len(seg)
            while k < len(ref) and ref[k] == base:
                run += 1; k += 1
            if run >= 5:
                hp = True
        for p in (1, 2, 3):                             # microssatélite curto
            if len(seg) >= 2 * p and len(seg) % p == 0 and seg == seg[:p] * (len(seg) // p):
                hp = True
    return '+'.join(segs), hp
