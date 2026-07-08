"""Tests for funvip.src.cluster._majority_top_group (the group tie-break)."""
import pandas as pd

from funvip.src.cluster import _majority_top_group


def _df(rows):
    return pd.DataFrame(rows, columns=["bitscore", "subject_group"]).sort_values(
        "bitscore", ascending=False
    )


def test_unique_top_bitscore_group_wins():
    # A alone holds the top bitscore -> A, regardless of lower tiers.
    assert _majority_top_group(_df([(250, "A"), (240, "B"), (240, "B")])) == "A"


def test_plurality_at_top_level_wins():
    assert _majority_top_group(_df([(250, "A"), (250, "A"), (250, "B")])) == "A"


def test_tie_folds_in_next_bitscore_level():
    # A and B tie at 250; fold in 240 where A has more -> A.
    d = _df([(250, "A"), (250, "B"), (240, "A"), (240, "A")])
    assert _majority_top_group(d) == "A"


def test_full_tie_returns_highest_bitscore_hit():
    # Perfect tie at every level -> the highest-bitscore hit (deterministic).
    d = _df([(250, "A"), (250, "B"), (240, "A"), (240, "B")])
    assert _majority_top_group(d) in {"A", "B"}
    # first row after the descending sort is the tiebreak
    assert _majority_top_group(d) == d.iloc[0]["subject_group"]
