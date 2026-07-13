"""Tests for funvip.src.tool.get_genus_species and its genus-list handling."""
import os

import pytest

from funvip.src import tool
from funvip.src.tool import get_genus_species

GENUS_FILE = os.path.join(os.path.dirname(tool.__file__), "..", "data", "genus_line.txt")


@pytest.fixture()
def genus_list():
    return tuple(open(GENUS_FILE).read().splitlines())


def test_uninitialized_raises_valueerror():
    # Regression: an uninitialized genus_file gave a cryptic NameError, not this.
    tool._genus_list_cache.clear()
    tool.__dict__.pop("genus_file", None)
    with pytest.raises(ValueError):
        get_genus_species("Aspergillus terreus strain")


@pytest.mark.parametrize(
    "text,expected",
    [
        ("Aspergillus terreus strain ABC", ("Aspergillus", "terreus")),
        ("Penicillium citrinum internal transcribed spacer", ("Penicillium", "citrinum")),
        ("Talaromyces marneffei sp. 3", ("Talaromyces", "marneffei")),
        ("Fusarium sp.", ("Fusarium", "sp.")),
        ("sequence with no known genus here", ("", "")),
    ],
)
def test_get_genus_species_with_explicit_list(genus_list, text, expected):
    # genus_list is a tuple (get_genus_species is @lru_cache'd -> args must hash).
    assert get_genus_species(text, genus_list=genus_list) == expected


def test_cached_global_matches_explicit_list(genus_list):
    class _Path:
        genusdb = GENUS_FILE

    tool._genus_list_cache.clear()
    tool.initialize_path(_Path())
    text = "Cladosporium cladosporioides voucher X"
    assert get_genus_species(text) == get_genus_species(text, genus_list=genus_list)
    assert GENUS_FILE in tool._genus_list_cache  # read once, cached
