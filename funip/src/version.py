# Version management module
import sys
import subprocess
from importlib.metadata import version

# import


"""
version = {
    "FunVIP": "0.3.19.0.1.3",
    "BLASTn": "",
    "MMseqs2": "",
    "MAFFT": "",
    "TrimAl": "",
    "Gblocks": "0.91b",
    "FastTree": "",
    "IQTREE2": "",
    "RAxML": "",
}
"""


class Version:
    def __init__(self, opt, path):
        self.FunVIP = ""
        self.GenMine = ""
        self.BLASTn = ""
        self.MMseqs2 = ""
        self.MAFFT = ""
        self.trimAl = ""
        self.Gblocks = "0.91b"
        self.Modeltest_NG = ""
        self.FastTree = ""
        self.IQTREE2 = ""
        self.RAxML = ""

        # For windows platform
        if sys.platform == "win32":
            ### FunVIP
            self.FunVIP = version("FunIP")

            ### GenMine
            self.GenMine = version("GenMine")

            ### BLASTn
            CMD = [f"{path.sys_path}/external/BLAST_Windows/bin/blastn.exe", "-version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            ## Format
            # blastn: blastn: 2.12.0+
            #   Package: blast 2.12.0, build Jun  4 2021 03:25:07
            self.BLASTn = stdout_str.split("\n")[0].split(" ")[1].strip()
            # print("BLASTn", self.BLASTn)

            ### MMSeqs2
            CMD = [
                f"{path.sys_path}/external/mmseqs_Windows/mmseqs.bat",
                "-h",
            ]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.MMseqs2 = stdout_str.split("Version: ")[1].split("\n")[0].strip()
            # print("MMseqs2", self.MMseqs2)

            ### MAFFT
            CMD = [
                f"{path.sys_path}/external/MAFFT_Windows/mafft-win/mafft.bat",
                "--version",
            ]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            # MAFFT, output is on the stderr
            stderr_str = stderr.decode("utf-8")
            self.MAFFT = stderr_str.split("\n")[-2].split(" ")[0].strip()
            # print("MAFFT", self.MAFFT)

            ### TrimAl
            CMD = [
                f"{path.sys_path}/external/trimal.v1.4/trimAl/bin/trimal.exe",
                "--version",
            ]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.trimAl = stdout_str.split("\n")[1].split(" ")[1]
            # print("trimAl", self.trimAl)

            ### Modeltest-ng
            ## Not supported in Windows
            self.Modeltest_NG = "not supported"

            ### FastTree
            ## Also use stderr of FastTree
            CMD = [
                f"{path.sys_path}/external/FastTree_Windows/FastTree.exe",
                "-expert",
            ]

            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            stderr_str = stderr.decode("utf-8")
            self.FastTree = stderr_str.split(" ")[4]
            # print("FastTree", self.FastTree)

            ### IQTREE
            CMD = [
                f"{path.sys_path}/external/iqtree/bin/iqtree2.exe",
                "--version",
            ]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.IQTREE2 = stdout_str.split(" ")[3]
            # print("IQTREE2", self.IQTREE2)

            ### RAxML
            if opt.avx is True:
                CMD = [
                    f"{path.sys_path}/external/RAxML_Windows/raxmlHPC-PTHREADS-AVX2.exe",
                    "-v",
                ]
            else:
                CMD = [
                    f"{path.sys_path}/external/RAxML_Windows/raxmlHPC-PTHREADS-SSE3.exe",
                    "-v",
                ]

            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.RAxML = stdout_str.split("\n")[2].split(" ")[4]
            # print("RAxML", self.RAxML)

        # For apple silicon platform
        elif sys.platform == "darwin":
            ### BLASTn
            CMD = ["blastn", "-version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            ## Format
            # blastn: blastn: 2.12.0+
            #   Package: blast 2.12.0, build Jun  4 2021 03:25:07
            self.BLASTn = stdout_str.split("\n")[0].split(" ")[1].strip()
            # print("BLASTn", self.BLASTn)

            ### MMSeqs2
            CMD = ["mmseqs", "-h"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.MMseqs2 = stdout_str.split("Version: ")[1].split("\n")[0].strip()
            # print("MMseqs2", self.MMseqs2)

            ### MAFFT
            CMD = ["mafft", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            # MAFFT, output is on the stderr
            stderr_str = stderr.decode("utf-8")
            self.MAFFT = stderr_str.split("\n")[-2].split(" ")[0].strip()
            # print("MAFFT", self.MAFFT)

            ### TrimAl
            CMD = ["trimal", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.trimAl = stdout_str.split("\n")[1].split(" ")[1]
            # print("trimAl", self.trimAl)

            ### Modeltest-ng
            ## Not supported in Windows
            CMD = ["modeltest-ng", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout_str = stdout.decode("utf-8")

            print("Modeltest_NG - should finish this code")
            print(stdout)
            self.Modeltest_NG = "not supported"

            raise Exception

            ### FastTree
            ## Also use stderr of FastTree
            CMD = ["FastTree", "-expert"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stderr_str = stderr.decode("utf-8")
            self.FastTree = stderr_str.split(" ")[4]
            # print("FastTree", self.FastTree)

            ### IQTREE
            CMD = ["IQTree2", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.IQTREE2 = stdout_str.split(" ")[3]
            # print("IQTREE2", self.IQTREE2)

            ### RAxML
            if opt.avx is True:
                CMD = ["raxmlHPC-PTHREADS-AVX2", "-v"]
            else:
                CMD = ["raxmlHPC-PTHREADS-SSE3", "-v"]

            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.RAxML = stdout_str.split("\n")[2].split(" ")[4]
            # print("RAxML", self.RAxML)

        # For linux platform
        else:
            ### BLASTn
            CMD = ["blastn", "-version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            ## Format
            # blastn: blastn: 2.12.0+
            #   Package: blast 2.12.0, build Jun  4 2021 03:25:07
            self.BLASTn = stdout_str.split("\n")[0].split(" ")[1].strip()
            # print("BLASTn", self.BLASTn)

            ### MMSeqs2
            CMD = ["mmseqs", "-h"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.MMseqs2 = stdout_str.split("Version: ")[1].split("\n")[0].strip()
            # print("MMseqs2", self.MMseqs2)

            ### MAFFT
            CMD = ["mafft", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            # MAFFT, output is on the stderr
            stderr_str = stderr.decode("utf-8")
            self.MAFFT = stderr_str.split("\n")[-2].split(" ")[0].strip()
            # print("MAFFT", self.MAFFT)

            ### TrimAl
            CMD = ["trimal", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.trimAl = stdout_str.split("\n")[1].split(" ")[1]
            # print("trimAl", self.trimAl)

            ### Modeltest-ng
            ## Not supported in Windows
            CMD = ["modeltest-ng", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout_str = stdout.decode("utf-8")

            print("Modeltest_NG - should finish this code")
            print(stdout)
            self.Modeltest_NG = "not supported"

            raise Exception

            ### FastTree
            ## Also use stderr of FastTree
            CMD = ["FastTree", "-expert"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stderr_str = stderr.decode("utf-8")
            self.FastTree = stderr_str.split(" ")[4]
            # print("FastTree", self.FastTree)

            ### IQTREE
            CMD = ["IQTree2", "--version"]
            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.IQTREE2 = stdout_str.split(" ")[3]
            # print("IQTREE2", self.IQTREE2)

            ### RAxML
            if opt.avx is True:
                CMD = ["raxmlHPC-PTHREADS-AVX2", "-v"]
            else:
                CMD = ["raxmlHPC-PTHREADS-SSE3", "-v"]

            result = subprocess.Popen(
                CMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = result.communicate()
            stdout_str = stdout.decode("utf-8")
            self.RAxML = stdout_str.split("\n")[2].split(" ")[4]
