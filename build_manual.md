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
conda install -c conda-forge grayskull
conda install -c conda-forge packaging
conda install conda-build
conda install conda-verify
conda install anaconda-client
conda install git
grayskull pypi FunID // Check if version is as you expected during this process
anaconda login
conda config --set anaconda_upload no
conda config --add channels cwseo
conda-build ./funid -c conda-forge // copy tar.bz2 location when this one ends
conda install anaconda-project
anaconda upload {Build file location} // tar.bz2 in conda-build log
conda remove -n FunID-condabuild
```
