Introduction
=================

.. _reference_to_overview:

What is scTriangulate?
------------------------
scTriangulate is a python3 package designed to harmonize ``conflicting annotations`` in single cell genomics studies. 
While this approach can be applied to any situation in which the same dataset has multiple conflicting annotations, common challenges includes:

1. Integrate results from the same or multiple unsupervised clustering algorithms (i.e. Leiden, Seurat, SnapATAC) using different ``resolutions``.

2. Integrate results from both unsupervised and supervised (i.e. cellHarmony, Seurat label transfer) ``clustering algorithms``.

3. Integrate results from different ``reference atlases``.

4. Integrate labels from ``multi-modality`` single cell datasets (CITE-Seq, Multiome, TEA-Seq, ASAP-Seq, etc.).

.. image:: ./_static/schema_chop.png
   :height: 180px
   :width: 800px
   :align: center
   :target: target

scTriangulate enables the user to ``mix-and-match`` individual clustering results by leveraging customizable 
statistical measures of single cell cluster stability in combination with a cooperative game theory approach (`Shapley Value <https://en.wikipedia.org/wiki/Shapley_value>`_) 
to obtain an single optimal solution.

.. note::
    A typical scRNA-Seq dataset (10k cells) with four provided annotation-sets can run in ~10 minutes on a laptop. For larger datasets (100k cells) or multiome 
    (GEX + ATAC) with > 100k features (gene + peak), it is recommended to run the software in a high-performance compute environment.

Inputs and Outputs
---------------------
scTriangulate has no limitations on the modalities you are working with, so far we have provided supports for ``RNA``, ``ADT``, ``ATAC``, ``Splicing``, ``DNA mutation``,
All you need to do is to construct an ``anndata`` that has the features expression values associated with each modality (i.e. gene for RNA, surface 
protein for ADT, peak or bins or kmers or motifs for ATAC, PSI values for splicing, a binary matrix indicating whether a cell has mutation or not for mutation data, etc).
Furthermore, we are planning to support TCR/BCR and spatial data in the future, the key is still how to select informative features from these new modalities.

scTriangulate works seemlessly with the popular `scanpy <https://scanpy.readthedocs.io/en/stable/>`_ package. In addtion, we offer 
a myriad of preprocessing convenient functions to ease file conversion. Currently we accept following formats:

    * **Anndata** (.h5 & .h5ad), the annotations are the columns in adata.obs
    * **mtx**, and the annotation information should be supplied as an addtional txt file (see below example and :ref:`reference_to_add_annotation`)
    * **dense matrix**, txt expression matrix, and the annotations should be supplied as an addtional txt file (see below example and :ref:`reference_to_add_annotation`).

    .. csv-table:: annotation txt file
        :file: ./_static/annotation_txt.csv
        :widths: 10,10
        :header-rows: 1

Optionally, users can supply their own UMAP embeddings, Please refer to :ref:`reference_to_add_umap` function for details.

All of the intermediate outputs and final clustering results, plus interactive visualization, will be automatically named and saved to the user-defined
directory. Each function provides a `save` argument, which allows the users to modify this default behaviour. 

With that, feel free to jump to the :ref:`tutorials` to get a sense about how to run the program (super easy)!

