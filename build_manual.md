# test build
```
pip install ./
```

# pypi build
In new conda environment
```
conda create -n FunID-build python=3.10
conda activate FunID-build
pip install twine
python setup.py bdist_wheel --universal
python setup.py sdist
twine upload dist/FunID-{YOUR_VERSION}* 	// use current build number
conda deactivate
conda env remove -n FunID-build
```


# conda build
In new conda environment
```
conda create -n FunID-condabuild python=3.9
conda activate FunID-condabuild
conda install --yes -c conda-forge grayskull packaging
conda install --yes conda-build conda-verify anaconda-client git urllib3
grayskull pypi FunID // Check if version is as you expected during this step
anaconda login
conda config --set anaconda_upload no
conda config --add channels cwseo
conda-build ./funid -c conda-forge // copy tar.bz2 location when this one ends
conda install anaconda-project --yes
anaconda upload {Build file location} // tar.bz2 in conda-build log
conda deactivate
conda env remove -n FunID-condabuild
```
