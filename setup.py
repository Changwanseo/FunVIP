# FunVIP/setup.py
from setuptools import setup, find_packages
from sys import platform
import subprocess

name = "FunVIP"
__version__ = "0.3.19.0.1.13"


# Define common setup options
setup_options = dict(
    name=name,
    version=__version__,
    description="Fungal Validation & Identification Pipeline",
    author="Changwan Seo",
    author_email="wan101010@snu.ac.kr",
    url="https://github.com/Changwanseo/FunVIP",
    python_requires="<3.13, >3.8",
    packages=[
        "funvip",
        "funvip.src",
        "funvip.data",
        "funvip.test_dataset",
        "funvip.preset",
    ],
    install_requires=[
        "biopython==1.78",
        "ete3==3.1.3",
        "Cython",
        "dendropy",
        "GenMine>=1.0.13",
        "lxml",
        "matplotlib",
        "numpy",
        "openpyxl==3.0.9",
        "pandas==1.4.2",
        "plotly==5.9.0",
        "psutil",
        "pyyaml",
        "sip>=4.19.4",
        "scikit-learn",
        "scipy",
        "tabulate",
        "unidecode==1.2.0",
        "xlrd==2.0.1",
        "xlsxwriter",
        "xmltodict==0.12.0",
    ],
    zip_safe=False,
    entry_points={
        "console_scripts": ["FunID = funvip.main:main", "FunVIP = funvip.main:main"]
    },
    package_dir={"funvip": "funvip"},
    package_data={
        "funvip.data": ["*.xlsx", "*.txt"],
        "funvip.test_dataset": ["**"],
        "funvip.preset": ["**"],
    },
    include_package_data=True,
    license="GPL3",
)

# Modify setup options based on platform
if platform == "darwin":
    pass
elif platform == "win32":
    setup_options["packages"].append("funvip.external")
    setup_options["install_requires"].append("PyQt5>=5.9.2")
    setup_options["package_data"]["funvip.external"] = ["**"]
else:  # linux maybe
    setup_options["install_requires"].append("PyQt5>=5.9.2")


def run_tests():
    result = subprocess.run(["FunVIP", "-h"], capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)


# Change setup options by platform
if platform == "win32":
    setup_options["packages"].append("funvip.external")
    setup_options["package_data"]["funvip.external"] = ["**"]
else:
    pass

print(setup_options)


# Run setup
setup(**setup_options, cmdclass={"test": run_tests})
