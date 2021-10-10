Introduction
=================

.. _reference_to_overview:

What is scTriangulate?
------------------------
scTriangulate is a python package designed to solve ``"multi-annotation"`` problem. When analyzing single cell omics data, you may encounter situations
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
biologically meaningful metrics to assess cluster goodness, and `Shapley Value <https://en.wikipedia.org/wiki/Shapley_value>`_ from cooperative theory 
to attain a single stable solution.

.. note::
    A typical scRNA-Seq dataset (10k cells) with four provided annotation-sets can run in ~10 minutes in a laptop. For larger datasets (100k cells) or multiome 
    (GEX + ATAC) with > 100k features (gene + peak), it is recommended to run the program in the high-performance compute environment.

Inputs and Outputs
---------------------
scTriangulate is designed for h5ad file, it works seemlessly with popular scanpy packages if you are familiar with it. In addtion to that, we offer 
a myriad of preprocessing convenient functions to ease the file conversion process, currently we accept following format:

    * **Anndata** (.h5 & .h5ad), the annotations are the columns in adata.obs
    * **mtx**, annotations information should be supplied as addtional txt file (see below example and :ref:`reference_to_add_annotation`)
    * **dense matrix**, txt expression matrix, annotations should be aupplied as addtional txt file (see below example and :ref:`reference_to_add_annotation`).

    .. csv-table:: annotation txt file
        :file: ./_static/annotation_txt.csv
        :widths: 10,10
        :header-rows: 1

Optionally, users can supply their own umap embeddings, Please refer to :ref:`reference_to_add_umap` function for the details.

All the intermediate outputs and final clustering results, plus interactive visualizations, will be automatically named and saved to the user-defined
repository. Each function provide `save` argument which allows the users to modify this default behaviour. 

