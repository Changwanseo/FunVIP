# Version management + external-tool preflight
import sys
import shutil
import logging
import subprocess
from importlib.metadata import version as _pkg_version

from funvip.src.exceptions import ConfigError


def _probe(candidates, args, parse):
    """Return a version string for the first available command among `candidates`,
    or '' if none is found or the probe/parse fails. Never raises."""
    for cmd in candidates:
        if shutil.which(cmd) is None:
            continue
        try:
            result = subprocess.run(
                [cmd, *args], capture_output=True, text=True, timeout=30
            )
        except Exception:
            continue
        try:
            parsed = parse(result.stdout or "", result.stderr or "")
        except Exception:
            parsed = ""
        if parsed:
            return parsed.strip()
    return ""


class Version:
    """Best-effort version stamp for the report. Probes only succeed for tools
    that are installed; anything missing or unparseable is left as '' rather than
    crashing the run (a run only needs the tools for its selected methods, which
    the preflight verifies separately)."""

    def __init__(self, opt, path):
        win = sys.platform == "win32"
        ext = f"{path.sys_path}/external"

        def cand(linux, win_path):
            return [win_path] if win else linux

        try:
            self.FunVIP = _pkg_version("FunVIP")
        except Exception:
            self.FunVIP = ""
        try:
            self.GenMine = _pkg_version("GenMine")
        except Exception:
            self.GenMine = ""

        self.BLASTn = _probe(
            cand(["blastn"], f"{ext}/BLAST_Windows/bin/blastn.exe"),
            ["-version"],
            lambda o, e: o.split("\n")[0].split(" ")[1] if o.strip() else "",
        )
        self.MMseqs2 = _probe(
            cand(["mmseqs"], f"{ext}/mmseqs_Windows/mmseqs.bat"),
            ["-h"],
            lambda o, e: o.split("Version: ")[1].split("\n")[0] if "Version: " in o else "",
        )
        self.MAFFT = _probe(
            cand(["mafft"], f"{ext}/MAFFT_Windows/mafft-win/mafft.bat"),
            ["--version"],
            lambda o, e: e.split("\n")[-2].split(" ")[0] if e.strip() else "",
        )
        self.trimAl = _probe(
            cand(["trimal"], f"{ext}/trimal.v1.4/trimAl/bin/trimal.exe"),
            ["--version"],
            lambda o, e: o.split("\n")[1].split(" ")[1] if o.count("\n") >= 1 else "",
        )
        self.Gblocks = "0.91b"
        self.Modeltest_NG = (
            "not supported"
            if win
            else _probe(
                ["modeltest-ng"],
                ["--version"],
                lambda o, e: o.split("ModelTest-NG ")[1].split(" ")[0]
                if "ModelTest-NG " in o
                else "",
            )
        )
        self.FastTree = _probe(
            cand(["FastTree", "fasttree"], f"{ext}/FastTree_Windows/FastTree.exe"),
            ["-expert"],
            lambda o, e: e.split(" ")[4] if len(e.split(" ")) > 4 else "",
        )
        self.IQTREE2 = _probe(
            cand(["iqtree", "iqtree2"], f"{ext}/iqtree/bin/iqtree2.exe"),
            ["--version"],
            lambda o, e: o.split(" ")[3] if len(o.split(" ")) > 3 else "",
        )
        if win:
            raxml_cand = [
                f"{ext}/RAxML_Windows/raxmlHPC-PTHREADS-AVX2.exe"
                if opt.avx
                else f"{ext}/RAxML_Windows/raxmlHPC-PTHREADS-SSE3.exe"
            ]
        elif opt.avx:
            raxml_cand = [
                "raxmlHPC-PTHREADS-AVX2",
                "raxmlHPC-PTHREADS-SSE3",
                "raxmlHPC",
            ]
        else:
            raxml_cand = ["raxmlHPC-PTHREADS-SSE3", "raxmlHPC"]
        self.RAxML = _probe(
            raxml_cand,
            ["-v"],
            lambda o, e: o.split("\n")[2].split(" ")[4]
            if o.count("\n") >= 2 and len(o.split("\n")[2].split(" ")) > 4
            else "",
        )


def _needed_tools(opt):
    """(label, [command candidates], required) for the tools the selected methods
    will actually invoke."""
    tools = []

    search = str(opt.method.search).lower()
    if search == "blast":
        tools += [
            ("BLAST+ (blastn)", ["blastn"], True),
            ("BLAST+ (makeblastdb)", ["makeblastdb"], True),
        ]
    elif search == "mmseqs":
        tools += [("MMseqs2", ["mmseqs"], True)]

    tools += [("MAFFT", ["mafft"], True)]

    if opt.method.tcs is True:
        tools += [("T-COFFEE (TCS)", ["t_coffee"], False)]

    trim = str(opt.method.trim).lower()
    if trim == "trimal":
        tools += [("trimAl", ["trimal"], True)]
    elif trim == "gblocks":
        tools += [("Gblocks", ["Gblocks"], True)]

    modeltest = str(opt.method.modeltest).lower()
    if modeltest in ("modeltest-ng", "modeltestng"):
        tools += [("modeltest-ng", ["modeltest-ng"], True)]
    elif modeltest == "iqtree":
        tools += [("IQ-TREE (ModelFinder)", ["iqtree", "iqtree2"], True)]

    tree = str(opt.method.tree).lower()
    if tree == "fasttree":
        tools += [("FastTree", ["FastTree", "fasttree"], True)]
    elif tree == "iqtree":
        tools += [("IQ-TREE", ["iqtree", "iqtree2"], True)]
    elif tree == "raxml":
        cands = (
            ["raxmlHPC-PTHREADS-AVX2", "raxmlHPC-PTHREADS-SSE3", "raxmlHPC"]
            if opt.avx
            else ["raxmlHPC-PTHREADS-SSE3", "raxmlHPC"]
        )
        tools += [("RAxML", cands, True)]

    # de-duplicate by label (e.g. IQ-TREE can be both modeltest and tree)
    seen, out = set(), []
    for label, cands, req in tools:
        if label not in seen:
            seen.add(label)
            out.append((label, cands, req))
    return out


def preflight(opt):
    """Verify the external tools the selected methods need are available; raise a
    clear ConfigError listing what is missing. Optional tools (TCS) only warn.
    On Windows the tools are bundled, so this is a no-op."""
    if sys.platform == "win32":
        logging.info("Windows platform: using the bundled external tools")
        return

    present, missing = [], []
    for label, candidates, required in _needed_tools(opt):
        found = next((c for c in candidates if shutil.which(c)), None)
        if found:
            present.append(f"{label} [{found}]")
        elif required:
            missing.append((label, candidates))
        else:
            logging.warning(
                f"Optional tool {label} not found on PATH "
                f"({', '.join(candidates)}); that step will be skipped"
            )

    if present:
        logging.info("External tools found: " + ", ".join(present))

    if missing:
        detail = "\n".join(
            f"  - {label}: none of [{', '.join(cands)}] found on PATH"
            for label, cands in missing
        )
        raise ConfigError(
            "Required external tools for the selected methods are not on PATH:\n"
            f"{detail}\n"
            "Install them (e.g. `conda install -c bioconda blast mmseqs2 mafft "
            "trimal fasttree iqtree raxml modeltest-ng t-coffee`) or choose "
            "different --search / --trim / --modeltest / --tree methods."
        )
