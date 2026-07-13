"""End-to-end integration test on the bundled `terrei` dataset.

Marked `integration` (deselected by default; run with `pytest -m integration`).
Requires the external tools (mmseqs/mafft/fasttree) on PATH and an e-mail for the
bundled GenBank accessions via FUNVIP_TEST_EMAIL.

Asserts the pipeline's STABLE invariants rather than an exact byte match, because
mmseqs/mafft/fasttree are nondeterministic across runs (a few SPECIES_ASSIGNED
cells vary); species assignment is checked for presence, not exact value.
"""
import csv
import os
import shutil
import subprocess

import pytest

pytestmark = pytest.mark.integration


@pytest.mark.skipif(
    shutil.which("mmseqs") is None or shutil.which("mafft") is None,
    reason="external tools (mmseqs/mafft) not on PATH",
)
def test_terrei_end_to_end(tmp_path):
    email = os.environ.get("FUNVIP_TEST_EMAIL")
    if not email:
        pytest.skip("set FUNVIP_TEST_EMAIL to run the terrei integration test")

    outdir = tmp_path / "out"
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    proc = subprocess.run(
        [
            "FunVIP", "--test", "terrei", "--email", email,
            "--thread", "4", "--memory", "8G",
            "--outdir", str(outdir), "--runname", "terrei",
        ],
        env=env, capture_output=True, text=True, timeout=1200,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]

    result = outdir / "terrei" / "terrei.result.csv"
    assert result.exists(), "no result.csv produced"

    rows = list(csv.DictReader(open(result)))
    queries = [r for r in rows if r["DATATYPE"].strip().lower() == "query"]
    assert len(queries) == 79
    assert all(r["GROUP_ASSIGNED"] == "Aspergillus" for r in queries)
    assert all((r["SPECIES_ASSIGNED"] or "").strip() for r in queries)
