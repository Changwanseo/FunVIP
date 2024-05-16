# FunIP/setup.py
from setuptools import setup, find_packages
from sys import platform
import subprocess

name = "FunVIP"
__version__ = "0.3.19.0.1.6"
# release = "0.3.19.0.0.1"


## Default setup options
# PyQt5 should be independently installed for m1 architecture
if platform == "darwin":
    setup_options = dict(
        name=name,
        version=__version__,
        description="Fungal Validation & Identification Pipeline",
        author="Changwan Seo",
        author_email="wan101010@snu.ac.kr",
        url="https://github.com/Changwanseo/FunIP",
        python_requires="<3.13, >3.8",
        packages=[
            "funip",
            "funip.src",
            "funip.data",
            "funip.test_dataset",
            "funip.preset",
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
            "PyQt5>=5.9.2",
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

else:
    # Windows: python 3.12, 3.11, 3.10, 3.9 pass
    # Linux: python 3.10 pass
    setup_options = dict(
        name=name,
        version=__version__,
        description="Fungal Validation & Identification Pipeline",
        author="Changwan Seo",
        author_email="wan101010@snu.ac.kr",
        url="https://github.com/Changwanseo/FunIP",
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
            "PyQt5",
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
