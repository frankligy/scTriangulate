# scTriangulate

scTriangulate is a Python package to mix-and-match conflicting clustering results in single cell analysis, and generate reconciled 
clustering solutions.

## Tutorials

Check our full [documentation.](https://sctriangulate.readthedocs.io/en/latest/get_started.html)

## Overview

![schema](./image/schema.png)

It can potentially be used in an array of settings:

1. Running same unsupervised clustering (i.e. Leiden) algorithm using different resolutions.

2. Running unsupervised clustering using different algorithms.

3. Running reference mapping tools using different reference atlases.

4. Clustering labels from matched multi-modalities (RNA, ADT, ATAC, etc)

![schuma_chop](./image/schema_chop.png)


## Citation

scTriangulate will be presented in [2021 CZI Single-Cell Biology Annual Meeting](https://docs.google.com/document/d/142W5qCsXpv9CyyvQhmu_Re-Wf6mN7zZrAzarMMI6Se4/edit). 

A preprint will come out soon.

## Reproducibility

All the scripts for reproducing the results are available in [reproduce folder](https://github.com/frankligy/scTriangulate/tree/main/reproduce), along 
with all the necessary input files and intermediate outputs which are avaiable in [Synapse storage](https://www.synapse.org/#!Synapse:syn26320337/files/).

## Contact

Guangyuan(Frank) Li

li2g2@mail.uc.edu

PhD student, Biomedical Informatics

Cincinnati Children’s Hospital Medical Center(CCHMC)

University of Cincinnati, College of Medicine
