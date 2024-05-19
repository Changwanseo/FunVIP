# test build
```
pip install ./
```

# pypi build
In new conda environment
```
conda create -n FunVIP-build python=3.10
conda activate FunVIP-build
pip install twine setuptools-ext build
python -m build
twine upload dist/FunVIP-{YOUR_VERSION}* --config-file {your pypirc file} 	// use current build number
conda deactivate
conda env remove -n FunVIP-build
```


# conda build
In new conda environment
```
conda create -n FunVIP-condabuild python=3.9
conda activate FunVIP-condabuild
conda install --yes -c conda-forge grayskull packaging
conda install --yes conda-build conda-verify anaconda-client git urllib3=1.26.15
grayskull pypi FunVIP // Check if version is as you expected during this step
anaconda login
conda config --set anaconda_upload no
conda config --add channels funid
conda-build ./funvip -c conda-forge // copy tar.bz2 location when this one ends
conda install anaconda-project --yes
anaconda upload {Build file location} // tar.bz2 in conda-build log
conda deactivate
conda env remove -n FunID-condabuild
```
