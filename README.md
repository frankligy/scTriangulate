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

1. Step1: Check the [what the tool does](https:github.com). (1 min)
2. step2: Understand there are three modules in scTriangulate package: (2 mins)
    * [ScTriangulate Class](https:github.com) (core functionalities)
    * [preprocessing module](https:github.com) (flexible file/gene format conversion plus normalization options)
    * [colors module](https:github.com) (Allow publication-quality figure generation)
3. step3: follow the two workflow example ([single modality](https:github.com), [multi-modal](https:github.com) (30 mins)
4. step4: [Optional] Check the [principl](https:github.com) part to understand the philosophy of developing the tool.


## Installation

scTriangulate requires python >= 3.7, conda virtual environment is highly recommended.

```bash
pip install sctriangulate
```

From source,

```bash
git clone https://github.com/frankligy/scTriangulate
cd ./scTriangulate
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
