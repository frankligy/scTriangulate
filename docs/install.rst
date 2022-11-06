Installation
===============

scTriangulate requires python >= 3.7. This software has been extensively tested using the conda virtual environment::

    pip install sctriangulate


From source code::

    # method 1: Using pip (recommened)
    pip install git+https://github.com/frankligy/scTriangulate.git

    # method 2: Using setuptools 
    git clone https://github.com/frankligy/scTriangulate
    cd ./scTriangulate
    conda create -n sctriangulate_env python=3.7
    conda activate sctriangulate_env
    pip install --upgrade setuptools==57.5.0   
    python setup.py install
    # make sure setuptools <58, I tested setuptools=57.5.0


.. note::

    The above approaches will take care of all dependencies for you, for the information, scTriangulate depends on

    * squidpy ==1.2.0
    * gseapy ==0.10.4
    * scrublet ==0.2.3
    * yattag
    * anytree
    * mygene ==3.2.2

A minitest is included::

    cd .
    pytest

    # or
    cd ./test
    pytest mini_test.py -v
