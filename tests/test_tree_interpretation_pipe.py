"""Tests for reusing a single gene tree interpretation for the concatenated tree."""
import os
from types import SimpleNamespace

import pytest

try:
    from funvip.src import tree_interpretation_pipe as tip
except Exception as e:
    pytest.skip(f"tree_interpretation_pipe unavailable: {e}", allow_module_level=True)

from funvip.src.reporter import Singlereport

ALIGNMENT = ">HS1HE\nACGT-ACGT\n>HS2HE\nACGTTACGT\n"
TREE = "(HS1HE:0.1,HS2HE:0.2);\n"


def _fi(n):
    return SimpleNamespace(hash=f"HS{n}HE")


def _dataset(qr, og):
    return SimpleNamespace(list_qr_FI=qr, list_og_FI=og, list_db_FI=[])


def _setup(tmp_path, genes=("its",), concat_qr=None, concat_alignment=ALIGNMENT):
    path = SimpleNamespace(
        out_tree=(tmp_path / "07_Tree").as_posix(),
        out_alignment=(tmp_path / "05_Alignment").as_posix(),
    )
    os.makedirs(f"{path.out_tree}/hash")
    os.makedirs(f"{path.out_alignment}/hash")
    opt = SimpleNamespace(runname="R")

    qr, og = [_fi(1)], [_fi(2)]
    dict_dataset = {"G": {g: _dataset(qr, og) for g in genes}}
    dict_dataset["G"]["concatenated"] = _dataset(
        qr if concat_qr is None else concat_qr, og
    )

    alignments = {g: ALIGNMENT for g in genes}
    alignments["concatenated"] = concat_alignment
    for g, text in alignments.items():
        with open(f"{path.out_tree}/hash/hash_R_G_{g}.nwk", "w") as f:
            f.write(TREE)
        with open(f"{path.out_alignment}/hash/R_hash_trimmed_G_{g}.fasta", "w") as f:
            f.write(text)

    width = len(ALIGNMENT.split("\n")[1])
    V = SimpleNamespace(
        dict_dataset=dict_dataset,
        partition={"G": {"len": {genes[0]: width}, "order": [genes[0]]}},
    )
    return V, path, opt


def test_single_gene_group_with_identical_inputs_is_reused(tmp_path):
    V, path, opt = _setup(tmp_path)
    assert tip._single_gene_twin(V, path, opt, "G") == "its"


def test_multi_gene_group_is_not_reused(tmp_path):
    V, path, opt = _setup(tmp_path, genes=("its", "tef"))
    assert tip._single_gene_twin(V, path, opt, "G") is None


def test_differing_query_list_is_not_reused(tmp_path):
    V, path, opt = _setup(tmp_path, concat_qr=[_fi(1), _fi(3)])
    assert tip._single_gene_twin(V, path, opt, "G") is None


def test_differing_alignment_is_not_reused(tmp_path):
    V, path, opt = _setup(tmp_path, concat_alignment=ALIGNMENT.replace("ACGT-", "ACGTA"))
    assert tip._single_gene_twin(V, path, opt, "G") is None


def test_partition_not_covering_the_alignment_is_not_reused(tmp_path):
    V, path, opt = _setup(tmp_path)
    V.partition["G"]["len"]["its"] -= 1
    assert tip._single_gene_twin(V, path, opt, "G") is None


def test_twin_interpretation_copies_outputs_and_renames_gene(tmp_path):
    out_tree = (tmp_path / "07_Tree").as_posix()
    os.makedirs(out_tree)
    path = SimpleNamespace(out_tree=out_tree)
    opt = SimpleNamespace(runname="R")
    files = {
        "hash_R_G_its_original.svg": "hash svg",
        "R_G_its_original.svg": "svg",
        "R_G_its.nwk": "interpreted tree",
        "R_G_concatenated.nwk": "tree before interpretation",
    }
    for name, text in files.items():
        with open(f"{out_tree}/{name}", "w") as f:
            f.write(text)
    source = SimpleNamespace(group="G", gene="its", tree_name="its tree")

    twin = tip._twin_interpretation(source, path, opt)

    def read(name):
        with open(f"{out_tree}/{name}") as f:
            return f.read()

    assert read("hash_R_G_concatenated_original.svg") == "hash svg"
    assert read("R_G_concatenated_original.svg") == "svg"
    assert read("R_G_concatenated_original.nwk") == "tree before interpretation"
    assert read("R_G_concatenated.nwk") == "interpreted tree"
    assert (twin.gene, source.gene) == ("concatenated", "its")
    assert twin.tree_name.endswith("hash/hash_R_G_concatenated.nwk")


def test_twin_visualization_copies_svg_and_relabels_reports(tmp_path):
    out_tree = (tmp_path / "07_Tree").as_posix()
    os.makedirs(out_tree)
    with open(f"{out_tree}/R_G_its.svg", "w") as f:
        f.write("<svg/>")
    report = Singlereport()
    report.hash, report.gene, report.species_assigned = "HS1HE", "its", "Genus species"
    source = SimpleNamespace(group="G", gene="its")

    twin_reports = tip._twin_visualization(
        source,
        [report],
        SimpleNamespace(out_tree=out_tree),
        SimpleNamespace(runname="R"),
    )

    with open(f"{out_tree}/R_G_concatenated.svg") as f:
        assert f.read() == "<svg/>"
    assert [r.gene for r in twin_reports] == ["concatenated"]
    assert report.gene == "its"
    assert vars(twin_reports[0]) == {**vars(report), "gene": "concatenated"}
