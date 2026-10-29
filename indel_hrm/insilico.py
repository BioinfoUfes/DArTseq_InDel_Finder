"""Exact-match in-silico PCR screen.

A unique tag in the genome does not guarantee that the designed primer PAIR amplifies a
single product. This module scans the supplied reference for inward-facing exact primer
matches and counts the predicted products PER HAPLOTYPE — consistent with the
haplotype-aware specificity used elsewhere (a homeolog amplifying once in each haplotype
is expected; two products in the SAME haplotype flags a non-specific pair).

Limitations (intentional): exact matches only. It does not model primer mismatches,
mismatch thermodynamics, amplification efficiency, or secondary structure. One predicted
product is encouraging, not proof of specificity.
"""
from collections import defaultdict

_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def _find_all(seq, pat, cap=300):
    pos = []
    start = 0
    while True:
        p = seq.find(pat, start)
        if p < 0:
            break
        pos.append(p)
        start = p + 1
        if len(pos) >= cap:
            break
    return pos


def _count_oriented(seq, p1, p2, min_len, max_len):
    """Products where p1 (left) and p2 (already reverse-complemented, right-binding) face
    inward on this strand, with total product length in [min_len, max_len]."""
    P1 = _find_all(seq, p1)
    if not P1:
        return 0
    P2 = _find_all(seq, p2)
    if not P2:
        return 0
    n = 0
    L2 = len(p2)
    for i in P1:
        for j in P2:
            if j >= i:
                length = j + L2 - i
                if min_len <= length <= max_len:
                    n += 1
    return n


def count_products(genome_index, fwd, rev, min_len, max_len):
    """Return dict haplotype -> predicted product count, scanning both strands."""
    rc_rev = revcomp(rev)
    rc_fwd = revcomp(fwd)
    per_hap = defaultdict(int)
    seen = set()
    for cid, (seq, hap) in genome_index.items():
        if id(seq) in seen:          # skip the version-less alias
            continue
        seen.add(id(seq))
        # plus strand: fwd ...... rc(rev)
        per_hap[hap] += _count_oriented(seq, fwd, rc_rev, min_len, max_len)
        # minus strand (equivalent to scanning the plus strand for rev ...... rc(fwd))
        per_hap[hap] += _count_oriented(seq, rev, rc_fwd, min_len, max_len)
    return dict(per_hap)


def screen_markers(usable, genome_index, cfg):
    """Attach in-silico PCR results to each usable marker's first primer pair.

    Sets r['insilico'] = {per_hap, total, specific}. 'specific' = at least one product and
    no more than one product in any single haplotype.
    """
    for r in usable:
        pair = r["primers"]["pairs"][0]
        # size window from THIS pair's designed amplicon (the reference may carry either
        # allele), with tolerance — not from the global config range, which the relaxed
        # primer3 retry may have exceeded.
        la, ls = len(pair["amplicon_long"]), len(pair["amplicon_short"])
        lo = max(10, min(la, ls) - 20)
        hi = max(la, ls) + 20
        per_hap = count_products(genome_index, pair["left_seq"], pair["right_seq"], lo, hi)
        total = sum(per_hap.values())
        specific = total >= 1 and (max(per_hap.values()) <= 1 if per_hap else False)
        r["insilico"] = dict(per_hap=per_hap, total=total, specific=specific)
    return usable
