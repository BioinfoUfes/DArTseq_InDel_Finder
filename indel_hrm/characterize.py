"""Characterize each InDel candidate: per-sample state, discriminating power, replicate
consistency, homopolymer context and strength of the evidence.

Length-marker states (dominant SilicoDArT):
  A  = only tag A present       B  = only tag B present
  AB = both present             N  = none (double null)     miss = missing data
"""
import itertools
from .indel_finder import gap_segment_and_homopolymer


def classify(va, vb):
    def present(v):
        if v in ('1', '2'):
            return True
        if v == '0':
            return False
        return None
    pa, pb = present(va), present(vb)
    if pa is None or pb is None:
        return 'miss'
    if pa and not pb:
        return 'A'
    if pb and not pa:
        return 'B'
    if pa and pb:
        return 'AB'
    return 'N'


def characterize(cand, markers, cult):
    allpairs = list(itertools.combinations(range(len(cult)), 2))
    out = []
    for cd in cand:
        ma, mb = markers[cd['ka']], markers[cd['kb']]
        seg, hp = gap_segment_and_homopolymer(ma['seq'], mb['seq'], cd['gaps'])
        labels = {c: classify(ma['call'].get(c, '-'), mb['call'].get(c, '-')) for c in cult}
        incons = {c: bool(ma['incons'].get(c) or mb['incons'].get(c)) for c in cult}
        clean = {c: labels[c] for c in cult if labels[c] in ('A', 'B', 'AB')}
        allst = {c: labels[c] for c in cult if labels[c] != 'miss'}

        def nsep(stmap):
            n = 0
            for i, j in allpairs:
                ci, cj = cult[i], cult[j]
                if ci in stmap and cj in stmap and stmap[ci] != stmap[cj]:
                    n += 1
            return n

        out.append(dict(
            chrom=ma['chrom'], pos=ma['pos'], idA=ma['id'], idB=mb['id'],
            src=ma['src'] + '-' + mb['src'], indel=cd['indel'], seg=seg, hp=hp,
            ident=cd['ident'], exact_anchor=cd.get('exact_anchor', 0),
            labels=labels, incons=incons,
            sep_clean=nsep(clean), sep_all=nsep(allst),
            n_clean=len(clean), n_all=len(allst),
            incons_n=sum(incons.values()),
            crA=ma['cr'], crB=mb['cr'], repA=ma['rep'], repB=mb['rep'],
            seqA=ma['seq'], seqB=mb['seq']))
    return out, allpairs
