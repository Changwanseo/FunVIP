"""Tests for input-table validation in funvip.src.validate_input.input_table."""
import csv
import os
from types import SimpleNamespace

import pytest

from funvip.src import validate_input
from funvip.src.validate_input import input_table

SEQ = "ACGT" * 30


def _path(tmp_path):
    for d in ("db", "query", "GenMine"):
        (tmp_path / d).mkdir(exist_ok=True)
    return SimpleNamespace(
        genusdb=os.path.join(
            os.path.dirname(validate_input.__file__), "..", "data", "genus_line.txt"
        ),
        out_db=str(tmp_path / "db"),
        out_query=str(tmp_path / "query"),
        GenMine=str(tmp_path / "GenMine"),
    )


def _opt(email=""):
    return SimpleNamespace(level="genus", gene=["its"], email=email, tableformat="csv")


def _table(tmp_path, name, rows):
    p = tmp_path / name
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, quoting=csv.QUOTE_ALL)
        w.writerow(["id", "genus", "species", "its"])
        w.writerows(rows)
    return p.as_posix()


def _run(tmp_path, tables, datatype="db", email="", shared=None):
    kwargs = {}
    if shared is not None:
        kwargs["id_origin"] = shared[0]
    funinfo_dict, genmine_flag, warnings, errors = input_table(
        funinfo_dict={} if shared is None else shared[1],
        path=_path(tmp_path),
        opt=_opt(email),
        table_list=tables,
        datatype=datatype,
        **kwargs,
    )
    return funinfo_dict, genmine_flag, warnings, [e for e in errors if e]


@pytest.fixture
def no_genmine(monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("GenMine must not be called")

    monkeypatch.setattr(validate_input.subprocess, "call", fail)


def test_single_token_accessions_are_sent_to_genmine(tmp_path, monkeypatch):
    # Formats the strict patterns do not cover (2 letters + 8 digits, RefSeq with
    # 7 digits) must keep reaching GenMine as before.
    monkeypatch.setattr(validate_input, "sleep", lambda s: None)
    monkeypatch.setattr(validate_input.subprocess, "call", lambda *a, **k: 1)
    table = _table(
        tmp_path,
        "db.csv",
        [
            ("A1", "Aspergillus", "terreus", "MN123456"),
            ("A2", "Aspergillus", "terreus", "PQ12345678"),
            ("A3", "Aspergillus", "terreus", "NR_1234567.1"),
            ("A4", "Aspergillus", "terreus", " MN654321.1 "),
            ("A5", "Aspergillus", "terreus", SEQ),
        ],
    )
    with pytest.raises(Exception):
        _run(tmp_path, [table], email="test@example.org")

    with open(tmp_path / "GenMine" / "Accessions.txt") as f:
        sent = sorted(line.strip() for line in f)
    assert sent == ["MN123456", "MN654321.1", "NR_1234567.1", "PQ12345678"]


def test_free_text_and_accession_lists_stop_at_validation(tmp_path, no_genmine):
    table = _table(
        tmp_path,
        "db.csv",
        [
            ("A1", "Aspergillus", "terreus", SEQ),
            ("A2", "Aspergillus", "terreus", "ITS from MN123456 and MN123457"),
            ("A3", "Aspergillus", "terreus", "MN123456, MN123457"),
            ("A4", "Aspergillus", "terreus", "MN123456;MN123457"),
        ],
    )
    _, genmine_flag, _, errors = _run(tmp_path, [table])

    assert genmine_flag == 0
    rejected = [e for e in errors if "neither a single GenBank accession" in e]
    assert len(rejected) == 3
    assert any("its of line 1 " in e and "ITS from MN123456" in e for e in rejected)


def test_fasta_cell_with_accession_header_is_read_as_sequence(tmp_path, no_genmine):
    table = _table(
        tmp_path,
        "db.csv",
        [("A1", "Aspergillus", "terreus", f">MN123456.1 Aspergillus terreus\n{SEQ}")],
    )
    funinfo_dict, genmine_flag, _, errors = _run(tmp_path, [table])

    assert errors == []
    assert genmine_flag == 0
    assert funinfo_dict["A1"].seq["its"] == SEQ


@pytest.mark.parametrize(
    "first, second",
    [("LÖ21-04", "LO21-04"), ("GB‐0065942", "GB-0065942")],
)
def test_ids_colliding_after_unicode_replacement_stop_at_validation(
    tmp_path, no_genmine, first, second
):
    table = _table(
        tmp_path,
        "db.csv",
        [
            (first, "Aspergillus", "terreus", SEQ),
            (second, "Aspergillus", "terreus", "ACGT" * 25),
        ],
    )
    _, _, _, errors = _run(tmp_path, [table])

    collisions = [e for e in errors if "colliding with id" in e]
    assert len(collisions) == 1
    assert repr(first) in collisions[0] and repr(second) in collisions[0]


def test_ids_colliding_across_db_and_query_tables(tmp_path, no_genmine):
    shared = ({}, {})
    db = _table(tmp_path, "db.csv", [("LO21-04", "Aspergillus", "terreus", SEQ)])
    query = _table(tmp_path, "query.csv", [("LÖ21-04", "", "", SEQ)])

    _, _, _, db_errors = _run(tmp_path, [db], shared=shared)
    _, _, _, query_errors = _run(tmp_path, [query], datatype="query", shared=shared)

    assert db_errors == []
    assert any("colliding with id 'LO21-04' in table" in e for e in query_errors)


def test_identical_id_in_two_tables_is_still_merged(tmp_path, no_genmine):
    first = _table(tmp_path, "db1.csv", [("A1", "Aspergillus", "terreus", SEQ)])
    second = _table(tmp_path, "db2.csv", [("A1", "Aspergillus", "terreus", SEQ)])

    funinfo_dict, _, warnings, errors = _run(tmp_path, [first, second])

    assert errors == []
    assert list(funinfo_dict) == ["A1"]
    assert any("Duplicate id A1" in w for w in warnings)


def test_ids_without_letters_or_digits_are_rejected(tmp_path, no_genmine):
    table = _table(
        tmp_path,
        "db.csv",
        [
            ("A1", "Aspergillus", "terreus", SEQ),
            ("--", "Aspergillus", "terreus", SEQ),
            ("–", "Aspergillus", "terreus", SEQ),
            ("...", "Aspergillus", "terreus", SEQ),
        ],
    )
    _, _, _, errors = _run(tmp_path, [table])

    assert any("Empty id found, line [1, 2, 3]" in e for e in errors)
