"""Tests for funvip.src.reporter.Report.update_statistics (per-group summary)."""
import pandas as pd

from funvip.src.reporter import Report


def test_update_statistics_maps_status_per_group():
    r = Report()
    r.query_result = pd.DataFrame(
        {
            "DATATYPE": ["query"] * 6,
            "GROUP_ASSIGNED": ["A", "A", "A", "B", "B", "A"],
            "STATUS": ["match", "assigned", "conflict", "new species", "failed", "match"],
        }
    )
    r.update_statistics()
    s = r.statistics

    a = s["GROUP"].index("A")
    assert s["IDENTIFIED"][a] == 3  # match + assigned + match
    assert s["MISIDENTIFIED"][a] == 1  # conflict
    assert s["TOTAL"][a] == 4

    b = s["GROUP"].index("B")
    assert s["NEW SPECIES CANDIDATE"][b] == 1
    assert s["UNDETERMINED"][b] == 1
    assert s["TOTAL"][b] == 2

    t = s["GROUP"].index("TOTAL")
    assert s["TOTAL"][t] == 6
    assert s["IDENTIFIED"][t] == 3
    assert s["MISIDENTIFIED"][t] == 1


def test_update_statistics_empty_is_safe():
    r = Report()
    r.query_result = None
    r.update_statistics()  # must not raise
    assert r.statistics["GROUP"] == []
