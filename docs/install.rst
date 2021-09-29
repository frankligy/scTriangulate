Installation
===============

scTriangulate requires python >= 3.7, conda virtual environment is highly recommended::

    pip install sctriangulate


From source::

    git clone https://github.com/frankligy/scTriangulate
    cd ./scTriangulate
    python setup.py install
    # make sure setuptools <58, I tested setuptools=57.5.0


A minitest is included::

    cd ./test
    python -W ignore mini_test.py
    # should pass all tests
