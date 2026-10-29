"""Tests for InDel detection and alignment metrics."""
from indel_hrm.config import Config
from indel_hrm.indel_finder import (_indel_check_bio, _indel_check_difflib,
                                     find_indels, gap_segment_and_homopolymer)

FL = "ACGTACGTACGTACGTACGTAC"   # 22 bp left flank
FR = "GTCAGTCAGTCAGTCAGTCAGT"   # 22 bp right flank
INDEL = "TTGGAACCTT"            # 10 bp


def test_clean_indel_both_aligners():
    a = FL + INDEL + FR
    b = FL + FR
    for fn in (_indel_check_bio, _indel_check_difflib):
        match, mism, internal, gaps, ident, anchor = fn(a, b)
        assert internal == len(INDEL)
        assert len(gaps) == 1
        assert gaps[0][0] == "delete"
        assert ident == 1.0
        assert anchor >= 15          # exact flanks on both sides


def test_insert_direction():
    a = FL + FR
    b = FL + INDEL + FR
    _, _, internal, gaps, _, _ = _indel_check_bio(a, b)
    assert internal == len(INDEL)
    assert gaps[0][0] == "insert"


def test_terminal_gap_ignored():
    a = "GGGGGGG" + FL + FR
    b = FL + FR
    _, _, internal, gaps, _, _ = _indel_check_bio(a, b)
    assert internal == 0
    assert gaps == []


def test_gap_segment_recovered():
    a = FL + INDEL + FR
    b = FL + FR
    gaps = _indel_check_bio(a, b)[3]
    seg, hp = gap_segment_and_homopolymer(a, b, gaps)
    assert seg == INDEL
    assert hp is False


def test_homopolymer_flagged():
    poly = "AAAAAAAA"
    a = FL + poly + FR
    b = FL + FR
    gaps = _indel_check_bio(a, b)[3]
    _, hp = gap_segment_and_homopolymer(a, b, gaps)
    assert hp is True


def _marker(mid, pos, seq, cr=1.0, rep=1.0):
    return dict(id=mid, chrom="chr1", pos=pos, seq=seq, src="SIL",
                cr=cr, rep=rep, call={}, incons={})


def test_find_indels_single_candidate():
    cfg = Config()
    a = FL + INDEL + FR
    b = FL + FR
    markers = [_marker("tagA", 100, a), _marker("tagB", 100, b)]
    cand = find_indels(markers, cfg)
    assert len(cand) == 1
    assert cand[0]["indel"] == len(INDEL)
    assert cand[0]["exact_anchor"] >= 15


def test_find_indels_respects_coloc_window():
    cfg = Config()
    a = FL + INDEL + FR
    b = FL + FR
    markers = [_marker("tagA", 100, a), _marker("tagB", 100 + cfg.coloc_window + 5, b)]
    assert find_indels(markers, cfg) == []
