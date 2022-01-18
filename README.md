[![Documentation Status](https://readthedocs.org/projects/sctriangulate/badge/?version=latest)](https://sctriangulate.readthedocs.io/en/latest/?badge=latest)  [![Pypi](https://img.shields.io/pypi/v/sctriangulate?logo=PyPI)](https://pypi.org/project/sctriangulate/)  [![Downloads](https://pepy.tech/badge/sctriangulate)](https://pypi.org/project/sctriangulate/)  [![Stars](https://img.shields.io/github/stars/frankligy/scTriangulate)](https://github.com/frankligy/scTriangulate/stargazers)



# scTriangulate

scTriangulate is a Python package to mix-and-match conflicting clustering results in single cell analysis, and generate reconciled clustering solutions.

scTriangulate leverages cooperative game theory (Shapley Value) in conjunction with complimentary stability metrics (i.e. reassign score, TFIDF score and SCCAF score) to intelligently integrate clustering solutions from nearly unlimited sources. Applied to multimodal datasets, this approach highlights new cell populations and mechanisms underlying lineage diversity.

Please don't hesitate to reach out to me if you have any questions (contact down the page), I will be responsive.

## Tutorials and Installation

Check our [full documentation and step-by-step tutorials.](https://sctriangulate.readthedocs.io/en/latest/get_started.html)

## Overview

![schema](./image/schema.png)

It can be used in an array of settings:

1. Integrate results from the same or multiple unsupervised clustering algorithms (i.e. Leiden, Seurat, SnapATAC) using different resolutions.

2. Integrate results from both unsupervised and supervised (i.e. cellHarmony, Seurat label transfer) clustering algorithms.

3. Integrate results from different reference atlases.

4. Integrate labels from multi-modality single cell datasets (CITE-Seq, Multiome, TEA-Seq, ASAP-Seq, etc.).

![schuma_chop](./image/schema_chop.png)

## Citation

scTriangulate: Decision-level integration of multimodal single-cell data. BioRxiv. Oct 2021 (https://www.biorxiv.org/content/10.1101/2021.10.16.464640v1)

## Reproducibility

All scripts for reproducing the analyses in the preprint are available in the [reproduce folder](https://github.com/frankligy/scTriangulate/tree/main/reproduce), along 
with all the necessary input files and intermediate outputs which are avaiable in [Synapse storage](https://www.synapse.org/#!Synapse:syn26320337/files/).

## Contact

Guangyuan(Frank) Li

Email: li2g2@mail.uc.edu

PhD student, Biomedical Informatics

Cincinnati Children’s Hospital Medical Center(CCHMC)

University of Cincinnati, College of Medicine
