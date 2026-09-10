"""Tests for funvip.src.dataset.FunVIP_var.validate_alignments cleanup steps."""
import logging
import os

import pytest
from Bio import SeqIO

from funvip.src.dataset import Dataset, FunVIP_var
from funvip.src.validate_input import Funinfo


class _Method:
    tcs = False


class _Opt:
    runname = "test"
    terminate = False
    verbose = 2
    thread = 1
    method = _Method()


class _Path:
    def __init__(self, out_alignment):
        self.out_alignment = out_alignment


def _FI(_id, _hash, datatype, group="G"):
    FI = Funinfo()
    FI.original_id = _id
    FI.id = _id
    FI.hash = _hash
    FI.genus = "Genus"
    FI.ori_species = "species"
    FI.datatype = datatype
    FI.group = group
    FI.adjusted_group = group
    FI.seq = {"gene1": "ATGCAT"}
    return FI


def _write_alignment(out_alignment, records):
    os.makedirs(f"{out_alignment}/hash", exist_ok=True)
    os.makedirs(f"{out_alignment}/failed", exist_ok=True)
    with open(f"{out_alignment}/test_trimmed_G_gene1.fasta", "w") as fp:
        for _hash, seq in records:
            fp.write(f">{_hash}\n{seq}\n")


def _make_V(records):
    """Five FIs (3 db, 1 query, 1 outgroup) in one group/gene dataset."""
    db = [_FI("DB1", "HS1HE", "db"), _FI("DB2", "HS2HE", "db"), _FI("DB3", "HS3HE", "db")]
    qr = [_FI("QR_empty", "HS64628HE", "query")]
    og = [_FI("OG1", "HS4HE", "db")]

    V = FunVIP_var()
    V.list_FI = db + qr + og
    V.dict_hash_FI = {FI.hash: FI for FI in V.list_FI}
    V.dict_dataset = {
        "G": {"gene1": Dataset(gene="gene1", group="G", list_qr=qr, list_db=db, list_og=og)}
    }
    return V


@pytest.fixture
def aligned(tmp_path):
    """One sequence is all gaps after trimming; the rest fully overlap."""
    out_alignment = str(tmp_path)
    _write_alignment(
        out_alignment,
        [
            ("HS1HE", "ATGCAT"),
            ("HS2HE", "ATGCAT"),
            ("HS3HE", "ATGCAT"),
            ("HS64628HE", "------"),
            ("HS4HE", "ATGCAT"),
        ],
    )
    return _make_V([]), _Path(out_alignment), _Opt()


def test_empty_after_trim_sequence_is_removed_without_crashing(aligned, caplog):
    # Regression: the removal branch used to look up the never-populated
    # dict_hash_id and die with KeyError: 'HS64628HE'.
    V, path, opt = aligned
    with caplog.at_level(logging.WARNING):
        V.validate_alignments(path=path, opt=opt)

    dataset = V.dict_dataset["G"]["gene1"]
    assert [FI.hash for FI in dataset.list_qr_FI] == []
    assert [FI.hash for FI in dataset.list_db_FI] == ["HS1HE", "HS2HE", "HS3HE"]
    assert [FI.hash for FI in dataset.list_og_FI] == ["HS4HE"]


def test_removal_warning_reports_the_readable_id_not_the_hash(aligned, caplog):
    V, path, opt = aligned
    with caplog.at_level(logging.WARNING):
        V.validate_alignments(path=path, opt=opt)

    removals = [
        rec.message for rec in caplog.records if "removed from dataset" in rec.message
    ]
    assert removals == ["QR_empty removed from dataset G gene1. Please check the alignment and see the region is correct"]


def test_trimmed_alignment_is_rewritten_without_the_empty_sequence(aligned):
    V, path, opt = aligned
    V.validate_alignments(path=path, opt=opt)

    seq_list = list(
        SeqIO.parse(f"{path.out_alignment}/test_trimmed_G_gene1.fasta", "fasta")
    )
    assert [seq.id for seq in seq_list] == ["HS1HE", "HS2HE", "HS3HE", "HS4HE"]


def test_intact_alignment_keeps_every_sequence(tmp_path, caplog):
    out_alignment = str(tmp_path)
    _write_alignment(
        out_alignment,
        [
            ("HS1HE", "ATGCAT"),
            ("HS2HE", "ATGCAT"),
            ("HS3HE", "ATGCAT"),
            ("HS64628HE", "ATGCAT"),
            ("HS4HE", "ATGCAT"),
        ],
    )
    V = _make_V([])
    with caplog.at_level(logging.WARNING):
        V.validate_alignments(path=_Path(out_alignment), opt=_Opt())

    dataset = V.dict_dataset["G"]["gene1"]
    assert len(dataset.list_qr_FI) == 1
    assert len(dataset.list_db_FI) == 3
    assert len(dataset.list_og_FI) == 1
    assert not [
        rec.message for rec in caplog.records if "removed from dataset" in rec.message
    ]
