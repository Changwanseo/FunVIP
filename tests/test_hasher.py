"""Tests for funvip.src.hasher — hash encode/decode round-trips.

Guards the single-pass HS<n>HE decode against the previous N-alternative-regex
implementation across every newick/svg flag combination and the edge cases that
would break a naive rewrite (adjacent hashes, prefix-overlapping HS1HE/HS12HE,
tokens not in the dict, unicode/special-char values).
"""
import re

import pandas as pd
import pytest

from funvip.src import hasher
from funvip.src.hasher import newick_legal, svg_legal


def _reference_decode(hash_dict, content, newick=True, svg=False):
    """The original N-alternative-regex decode, kept as an oracle."""
    if newick and svg:
        hd = {re.escape(k): svg_legal(newick_legal(v)) for k, v in hash_dict.items()}
    elif newick:
        hd = {re.escape(k): newick_legal(v) for k, v in hash_dict.items()}
    elif svg:
        hd = {re.escape(k): svg_legal(v) for k, v in hash_dict.items()}
    else:
        hd = {re.escape(k): v for k, v in hash_dict.items()}
    pattern = re.compile("|".join(hd.keys()))
    return pattern.sub(lambda m: hd[re.escape(m.group(0))], content)


HASH_DICT = {
    "HS0HE": "AB123 Aspergillus terreus (strain A):var. B; 'note' & <x>",
    "HS1HE": "CD456 Penicillium citrinum",
    "HS12HE": "EF789 Talaromyces_marneffei, sp. 3",
    "HS2HE": "GH000 Fungus ésp. 中文",
    "HS10HE": "IJ111 Genus species",
}
CONTENT = (
    "(HS0HE:0.1,(HS1HE:0.2,HS12HE:0.3):0.05,HS2HE:0.4);\n"
    "adjacent:HS1HE HS0HEHS1HE end\n"
    "overlap:HS12HE vs HS1HE vs HS10HE\n"
    "not-in-dict:HS999HE and HS1HEZ trailing\n"
    "plain line with no hashes\n"
)


@pytest.mark.parametrize("newick,svg", [(True, True), (True, False), (False, True), (False, False)])
def test_decode_matches_reference(tmp_path, newick, svg):
    infile = tmp_path / "in.txt"
    infile.write_text(CONTENT, encoding="utf-8")
    out = tmp_path / "out.txt"
    hasher.decode(HASH_DICT, str(infile), str(out), newick=newick, svg=svg)
    assert out.read_text(encoding="utf-8") == _reference_decode(HASH_DICT, CONTENT, newick=newick, svg=svg)


def test_decode_leaves_unknown_hash(tmp_path):
    infile = tmp_path / "in.txt"
    infile.write_text("keep HS999HE unchanged\n", encoding="utf-8")
    out = tmp_path / "out.txt"
    hasher.decode({"HS0HE": "x"}, str(infile), str(out), newick=False)
    assert out.read_text(encoding="utf-8") == "keep HS999HE unchanged\n"


def test_decode_df_exact_cells_only():
    df = pd.DataFrame(
        {"qseqid": ["HS0HE", "HS12HE", "notahash"], "sseqid": ["HS1HE", "HS999HE", "HS2HE"]}
    )
    out = hasher.decode_df(HASH_DICT, df)
    assert list(out["qseqid"]) == [HASH_DICT["HS0HE"], HASH_DICT["HS12HE"], "notahash"]
    assert list(out["sseqid"]) == [HASH_DICT["HS1HE"], "HS999HE", HASH_DICT["HS2HE"]]


def test_decode_df_does_not_mutate_input():
    df = pd.DataFrame({"qseqid": ["HS0HE"]})
    hasher.decode_df(HASH_DICT, df)
    assert list(df["qseqid"]) == ["HS0HE"]
