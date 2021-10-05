# scTriangulate

scTriangulate is a Python package to mix-and-match conflicting clustering results in single cell analysis, and generate reconciled 
clustering solutions.

## Overview

![schema](./image/schema.png)

It can potentially be used in an array of settings:

1. Running same unsupervised clustering (i.e. Leiden) algorithm using different resolutions.

2. Running unsupervised clustering using different algorithms.

3. Running reference mapping tools using different reference atlases.

4. Clustering labels from matched multi-modalities (RNA, ADT, ATAC, etc)

![schuma_chop](./image/schema_chop.png)

## Get started

Check our documentation.


## Installation

scTriangulate requires python >= 3.7, conda virtual environment is highly recommended.

```bash
pip install sctriangulate
```

From source,

```bash
git clone https://github.com/frankligy/scTriangulate
cd ./scTriangulate
conda create -n sctriangulate_env python=3.7
conda activate sctriangulate_env
pip install --upgrade setuptools==57.5.0
python setup.py install
# make sure setuptools <58, I tested setuptools=57.5.0
```

A minitest is included,

```bash
cd ./test
python -W ignore mini_test.py
```

## Contact

Guangyuan(Frank) Li

li2g2@mail.uc.edu

PhD student, Biomedical Informatics

Cincinnati Children’s Hospital Medical Center(CCHMC)

University of Cincinnati, College of Medicine
