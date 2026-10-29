"""Tests for genotype classification and genotype-aware separation (incl. N handling)."""
import numpy as np

from indel_hrm.characterize import classify
from indel_hrm.genotype import pair_separation


def test_classify_states():
    assert classify("1", "0") == "A"
    assert classify("0", "1") == "B"
    assert classify("1", "1") == "AB"
    assert classify("0", "0") == "N"
    assert classify("-", "1") == "miss"
    assert classify("1", "-") == "miss"


def _marker():
    return dict(
        curves={"A": np.linspace(0, 1, 20),
                "B": np.linspace(1, 0, 20),
                "AB": np.full(20, 0.5)},
        labels={"x": "A", "y": "B", "z": "N", "w": "AB", "x2": "A", "m": "miss"},
    )


def test_A_vs_B_separates_strong():
    peak, rms, strong = pair_separation(_marker(), "x", "y")
    assert peak > 0
    assert strong is True


def test_same_genotype_no_separation():
    peak, rms, strong = pair_separation(_marker(), "x", "x2")
    assert peak == 0.0
    assert strong is True          # same state: truly identical curve, not an N problem


def test_N_not_strong_and_zero():
    # N has no predicted curve: must not be scored as a real separation
    peak, rms, strong = pair_separation(_marker(), "x", "z")
    assert peak == 0.0
    assert strong is False


def test_missing_not_strong():
    peak, rms, strong = pair_separation(_marker(), "x", "m")
    assert peak == 0.0
    assert strong is False


def test_AB_vs_A_separates():
    peak, rms, strong = pair_separation(_marker(), "x", "w")
    assert peak > 0
    assert strong is True
