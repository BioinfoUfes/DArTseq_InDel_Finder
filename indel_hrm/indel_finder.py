"""Detect InDel candidates from co-located DArTseq tags.

Because DArTseq-SNP only reports substitutions, InDels are inferred: we group tags that
DArT's BLAST placed at the same locus (same contig, positions <= coloc_window bp apart)
and align each pair; an INTERNAL (non-terminal) gap of size indel_min..indel_max, with
high-identity flanks, is an InDel candidate.

Alignment uses a proper DNA pairwise aligner with affine gap penalties (Biopython's
PairwiseAligner, free end gaps = semi-global) so that match/mismatch/gap scoring is
explicit and defensible. A cheap difflib pre-gate avoids aligning obvious non-pairs, and
a difflib fallback keeps the module working if Biopython is unavailable.
"""
from collections import defaultdict
from difflib import SequenceMatcher

_ALIGNER = None


def _get_aligner():
    """Lazily build a semi-global affine-gap DNA aligner (free terminal gaps)."""
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner
        a = PairwiseAligner()
        a.mode = "global"
        a.match_score = 1.0
        a.mismatch_score = -1.0
        a.open_gap_score = -2.0
        a.extend_gap_score = -0.5
        a.end_gap_score = 0.0          # free overhang at the ends (trimming)
        _ALIGNER = a
    return _ALIGNER


def _metrics_from_rows(rowA, rowB, seqA, seqB):
    """Walk the two gapped alignment rows and return the same contract as before:
    (match, mism, internal, gaps, ident, exact_anchor).

    gaps: list of (tag, size, i1, i2, j1, j2) in ORIGINAL coordinates, where
    tag='delete' => bases in seqA[i1:i2] (absent from B), tag='insert' => seqB[j1:j2].
    Terminal gap runs (leading/trailing overhang) are ignored. exact_anchor is the
    smallest exact (identical) matched run directly flanking an internal gap.
    """
    cols = len(rowA)
    nonrun = [i for i in range(cols) if rowA[i] != "-" and rowB[i] != "-"]
    if not nonrun:
        return 0, 0, 0, [], 0.0, 0
    lead, trail = nonrun[0], nonrun[-1]

    match = mism = internal = 0
    gaps = []                      # (tag, g, i1, i2, j1, j2, left_anchor, end_col)
    ia = ib = 0
    cur_run = 0                    # identical-match run length ending at the current point
    i = 0
    while i < cols:
        ca, cb = rowA[i], rowB[i]
        if ca != "-" and cb != "-":
            if ca == cb:
                match += 1
                cur_run += 1
            else:
                mism += 1
                cur_run = 0
            ia += 1
            ib += 1
            i += 1
            continue
        # gap run: group consecutive columns of the same gap type
        terminal = i < lead or i > trail
        la = cur_run                        # exact-match run ending just before this gap
        if ca == "-":                       # bases in B only -> insert
            j1 = ib
            while i < cols and rowA[i] == "-" and rowB[i] != "-":
                ib += 1
                i += 1
            g = ib - j1
            if not terminal:
                internal += g
                gaps.append(("insert", g, ia, ia, j1, ib, la, i))
        else:                               # bases in A only -> delete
            i1 = ia
            while i < cols and rowB[i] == "-" and rowA[i] != "-":
                ia += 1
                i += 1
            g = ia - i1
            if not terminal:
                internal += g
                gaps.append(("delete", g, i1, ia, ib, ib, la, i))
        cur_run = 0

    # exact flanking anchors around each internal gap
    exact_anchor = 10 ** 9
    for gp in gaps:
        left, end_col = gp[6], gp[7]
        # right exact run: count identical matched columns starting at end_col
        right = 0
        k = end_col
        while k < cols and rowA[k] != "-" and rowB[k] != "-" and rowA[k] == rowB[k]:
            right += 1
            k += 1
        exact_anchor = min(exact_anchor, left, right)
    exact_anchor = 0 if exact_anchor == 10 ** 9 else exact_anchor

    aligned = match + mism
    ident = match / aligned if aligned else 0.0
    # strip the helper fields from the gap tuples (keep the classic 6-field shape)
    gaps = [(t, g, i1, i2, j1, j2) for (t, g, i1, i2, j1, j2, _la, _ec) in gaps]
    return match, mism, internal, gaps, ident, exact_anchor


def _indel_check_bio(a, b):
    """Affine-gap alignment; returns (match, mism, internal, gaps, ident, exact_anchor)."""
    aln = _get_aligner().align(a, b)[0]
    rowA, rowB = str(aln[0]), str(aln[1])
    return _metrics_from_rows(rowA, rowB, a, b)


def _indel_check_difflib(a, b):
    """Fallback aligner (text diff). Same contract as _indel_check_bio."""
    sm = SequenceMatcher(None, a, b, autojunk=False)
    ops = sm.get_opcodes()
    match = mism = internal = 0
    gaps = []
    n = len(ops)
    exact_anchor = 10 ** 9
    for oi, (tag, i1, i2, j1, j2) in enumerate(ops):
        if tag == "equal":
            match += (i2 - i1)
        elif tag == "replace":
            mism += max(i2 - i1, j2 - j1)
        elif tag in ("delete", "insert"):
            g = (i2 - i1) + (j2 - j1)
            if oi not in (0, n - 1):
                internal += g
                gaps.append((tag, g, i1, i2, j1, j2))
                left = (ops[oi - 1][2] - ops[oi - 1][1]) if oi > 0 and ops[oi - 1][0] == "equal" else 0
                right = (ops[oi + 1][2] - ops[oi + 1][1]) if oi < n - 1 and ops[oi + 1][0] == "equal" else 0
                exact_anchor = min(exact_anchor, left, right)
    exact_anchor = 0 if exact_anchor == 10 ** 9 else exact_anchor
    aligned = match + mism
    ident = match / aligned if aligned else 0.0
    return match, mism, internal, gaps, ident, exact_anchor


def _quick_gate(a, b, cfg):
    """Cheap difflib screen: is there plausibly an internal indel in range? Avoids
    running the full aligner on obvious non-pairs."""
    sm = SequenceMatcher(None, a, b, autojunk=False)
    ops = sm.get_opcodes()
    n = len(ops)
    g = 0
    for oi, (tag, i1, i2, j1, j2) in enumerate(ops):
        if tag in ("delete", "insert") and oi not in (0, n - 1):
            g += (i2 - i1) + (j2 - j1)
    return cfg.indel_min <= g <= cfg.indel_max * 2 + 5


def _indel_check(a, b, cfg):
    if _quick_gate(a, b, cfg):
        try:
            return _indel_check_bio(a, b)
        except Exception:
            return _indel_check_difflib(a, b)
    return 0, 0, 0, [], 0.0, 0


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
                match, mism, internal, gaps, ident, anchor = _indel_check(ma['seq'], mb['seq'], cfg)
                if internal < cfg.indel_min or internal > cfg.indel_max:
                    continue
                if ident < cfg.min_ident or match < cfg.min_anchor:
                    continue
                if len(gaps) > cfg.max_gap_runs:
                    continue
                if cfg.min_exact_anchor and anchor < cfg.min_exact_anchor:
                    continue
                cand.append(dict(ka=ka, kb=kb, match=match, mism=mism,
                                 indel=internal, gaps=gaps, ident=round(ident, 3),
                                 exact_anchor=anchor))
    return cand


def gap_segment_and_homopolymer(seqA, seqB, gaps):
    """Return (InDel sequence, homopolymer/STR flag)."""
    hp = False
    segs = []
    for tag, g, i1, i2, j1, j2 in gaps:
        if tag == 'delete':
            seg = seqA[i1:i2]; pos = i1; ref = seqA
        else:
            seg = seqB[j1:j2]; pos = j1; ref = seqB
        segs.append(seg)
        if seg and len(set(seg)) == 1:                  # homopolymer
            base = seg[0]; run = len(seg); k = pos - 1
            while k >= 0 and ref[k] == base:
                run += 1; k -= 1
            k = pos + len(seg)
            while k < len(ref) and ref[k] == base:
                run += 1; k += 1
            if run >= 5:
                hp = True
        for p in (1, 2, 3):                             # short microsatellite
            if len(seg) >= 2 * p and len(seg) % p == 0 and seg == seg[:p] * (len(seg) // p):
                hp = True
    return '+'.join(segs), hp
