"""Tests for the exact-match in-silico PCR screen."""
from indel_hrm.insilico import revcomp, count_products


def test_revcomp():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAAA") == "TTTT"
    assert revcomp("ACCGGT") == "ACCGGT"


FWD = "ACGTACGTACGTACGTACGT"        # 20 bp
RC_REV = "TTTTGGGGCCCCAAAATTTT"     # 20 bp (right-binding site on the + strand)
REV = revcomp(RC_REV)


def test_single_product_counted():
    mid = "TTTTTTTTTTTTTTTTTTTT"    # 20 bp -> product length 60
    seq = "G" * 40 + FWD + mid + RC_REV + "C" * 40
    per_hap = count_products({"c1": (seq, "hapX")}, FWD, REV, 40, 130)
    assert per_hap.get("hapX", 0) == 1


def test_no_product_out_of_size_range():
    mid = "T" * 300                  # product too long
    seq = FWD + mid + RC_REV
    per_hap = count_products({"c1": (seq, "hapX")}, FWD, REV, 40, 130)
    assert per_hap.get("hapX", 0) == 0


def test_two_products_same_haplotype_flagged():
    mid = "T" * 20
    one = FWD + mid + RC_REV
    seq = "G" * 30 + one + "A" * 30 + one + "C" * 30
    per_hap = count_products({"c1": (seq, "hapX")}, FWD, REV, 40, 130)
    assert per_hap.get("hapX", 0) >= 2
