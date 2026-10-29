"""Tests for ranking and the minimal-panel greedy set cover."""
from indel_hrm.config import Config
from indel_hrm.select import rank, minimal_panel


def _cand(labels, indel=10, ident=0.95, sep_clean=None):
    import itertools
    if sep_clean is None:
        keys = list(labels)
        sep_clean = sum(1 for a, b in itertools.combinations(keys, 2)
                        if labels[a] in ("A", "B", "AB") and labels[b] in ("A", "B", "AB")
                        and labels[a] != labels[b])
    return dict(src="SIL-SIL", hp=False, indel=indel, incons_n=0, ident=ident,
                sep_clean=sep_clean, labels=labels)


def test_minimal_panel_covers_all_pairs():
    targets = ["a", "b", "c"]
    out = [
        _cand({"a": "A", "b": "B", "c": "B"}),   # separates (a,b),(a,c)
        _cand({"a": "A", "b": "A", "c": "B"}),   # separates (a,c),(b,c)
    ]
    cfg = Config()
    mp = minimal_panel(out, targets, cfg)
    assert mp["n_pairs"] == 3
    assert mp["resolvable"] == 3
    assert mp["resolved"] == 3
    assert mp["unresolved"] == []
    assert 1 <= len(mp["panel"]) <= 2


def test_unresolvable_pair_reported():
    targets = ["a", "b", "c"]
    # no marker ever separates b from c (always same state)
    out = [_cand({"a": "A", "b": "B", "c": "B"})]
    cfg = Config()
    mp = minimal_panel(out, targets, cfg)
    assert ("b", "c") in [tuple(p) for p in mp["unresolved"]]


def test_rank_prefers_higher_discrimination():
    low = _cand({"a": "A", "b": "A", "c": "A"})    # separates nothing -> sep_clean 0
    high = _cand({"a": "A", "b": "B", "c": "AB"})  # separates all 3 pairs -> sep_clean 3
    ranked = rank([low, high], Config())
    assert ranked[0] is high
    assert ranked[0]["sep_clean"] > ranked[-1]["sep_clean"]
