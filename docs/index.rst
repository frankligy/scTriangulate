.. scTriangulate documentation master file, created by
   sphinx-quickstart on Thu Jul 22 16:24:43 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



What is scTriangulate?
==========================

scTriangulate is a python package designed to solve ``"multi-annotation"`` problem. When analyzing single cell data, you may encounter situations
where you inevitably have multiple different sets of cell annotations:

1. Running same unsupervised clustering (i.e. Leiden) algorithm using different ``resolutions``.

2. Running unsupervised clustering using different ``algorithms``.

3. Running reference mapping tools using different ``reference atlases``.

4. Clustering labels from matched ``multi-modalities`` (RNA, ADT, ATAC, etc)

.. image:: ./_static/schema_chop.png
   :height: 180px
   :width: 800px
   :align: center
   :target: target

scTriangulate enable the usrs to ``mix-and-match`` the individual clustering results via leveraging customizable 
biologically meaningful metrics to assess cluster goodness, and ``Shapley Value`` from cooperative theory to attain a single stable solution.


Contents
===============

.. toctree::
   :maxdepth: 2

   install
   tutorial
   principle
   api
   contact



