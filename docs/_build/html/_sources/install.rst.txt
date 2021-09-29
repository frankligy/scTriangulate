Installation
===============

scTriangulate requires python >= 3.7, conda virtual environment is highly recommended::

    pip install sctriangulate


From source::

    git clone https://github.com/frankligy/scTriangulate
    cd ./scTriangulate
    python setup.py install
    # make sure setuptools <58, I tested setuptools=57.5.0


.. note::

    The above approaches will take care of all dependencies for you, for the information, scTriangulate depends on

    * 'scanpy ==1.7.2',
    * 'gseapy ==0.10.4',
    * 'scrublet ==0.2.3',
    * 'yattag',
    * 'anytree',
    * 'mygene ==3.2.2',
    * 'numpy <1.21'




A minitest is included::

    cd ./test
    python -W ignore mini_test.py
    # should pass all tests
