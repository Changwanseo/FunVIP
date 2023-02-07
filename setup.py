# FunID/setup.py
from setuptools import setup, find_packages
from sys import platform
import subprocess

name = "FunID"
__version__ = "0.3.1"
release = "0.3.1.0.0"

# Default setup options
setup_options = dict(
    name=name,
    version=__version__,
    description="Fungal Identification Pipeline",
    author="Changwan Seo",
    author_email="wan101010@snu.ac.kr",
    url="https://github.com/Changwanseo/FunID",
    python_requires="<3.11, >3.8",
    packages=["funid", "funid.src", "funid.data", "funid.test_dataset", "funid.preset"],
    install_requires=[
        "biopython==1.78",
        "ete3==3.1.2",
        "Cython",
        "dendropy",
        "GenMine",
        "lxml",
        "matplotlib==3.5.1",
        "numpy==1.22.3",
        "openpyxl==3.0.9",
        "pandas==1.4.2",
        "plotly==5.9.0",
        "PyQt5>=5.9.2",
        "pyyaml",
        "sip>=4.19.4",
        "scikit-learn==1.0.2",
        "tabulate",
        "unidecode==1.2.0",
        "xlrd==2.0.1",
        "xlsxwriter",
        "xmltodict==0.12.0",
    ],
    zip_safe=False,
    entry_points={"console_scripts": ["FunID = funid.main:main"]},
    package_dir={"funid": "funid"},
    package_data={
        "funid.data": ["*.xlsx", "*.txt"],
        "funid.test_dataset": ["**"],
        "funid.preset": ["**"],
    },
    include_package_data=True,
    license="GPL3",
)


def run_tests():
    result = subprocess.run(["FunID", "-h"], capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)


# Change setup options by platform
if platform == "win32":
    setup_options["packages"].append("funid.external")
    setup_options["package_data"]["funid.external"] = ["**"]
else:
    pass

print(setup_options)


# Run setup
setup(**setup_options, cmdclass={"test": run_tests})
