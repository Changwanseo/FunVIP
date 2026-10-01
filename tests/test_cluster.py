"""Tests for funvip.src.cluster: the group tie-break and outgroup task generation."""
from types import SimpleNamespace

import pandas as pd

from funvip.src.cluster import _majority_top_group, outgroup_append_opt_generator


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


def test_outgroup_tasks_group_the_search_table_once(monkeypatch):
    calls = []
    groupby = pd.DataFrame.groupby

    def counting_groupby(self, *args, **kwargs):
        calls.append(1)
        return groupby(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "groupby", counting_groupby)
    cSR = pd.DataFrame(
        {
            "qseqid": ["HS1HE", "HS2HE", "HS3HE", "HS4HE"],
            "bitscore": [300.0, 250.0, 200.0, 150.0],
            "query_group": ["A", "B", "A", "C"],
        }
    )
    V = SimpleNamespace(
        cSR=cSR,
        dict_dataset={
            "A": {"its": None, "concatenated": None},
            "B": {"its": None, "tef": None, "concatenated": None},
            "C": {"its": None},
            "D": {"its": None, "concatenated": None},
        },
    )

    tasks = outgroup_append_opt_generator(V, None, None)

    assert len(calls) == 1
    assert [(t[2], t[1]) for t in tasks] == [
        ("A", "its"),
        ("A", "concatenated"),
        ("B", "its"),
        ("B", "tef"),
        ("B", "concatenated"),
    ]
    assert list(tasks[0][0]["qseqid"]) == ["HS1HE", "HS3HE"]
    assert list(tasks[2][0]["qseqid"]) == ["HS2HE"]
