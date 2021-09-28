Tutorial
==========

Get started
------------------------

Hi, Welcome to scTriangulate tutorials, here are some suggestive steps to quickly get a hang of this tool:

1. Step1: Check the `what the tool does <https:github.com>`_. (1 min)
2. step2: Understand there are three modules in scTriangulate package: (2 mins)
    * `ScTriangulate Class <https:github.com>`_ (core functionalities)
    * `preprocessing module <https:github.com>`_ (flexible file/gene format conversion plus normalization options)
    * `colors module <https:github.com>`_ (Allow publication-quality figure generation)
3. step3: follow the two workflow example (`single modality <https:github.com>`_, `multi-modal <https:github.com>`_) (30 mins)
4. step4: [Optional] Check the `principle <https:github.com>`_ part to understand the philosophy of developing the tool.



Single Modality (scRNA) workflow
-----------------------------------

In this example, we are going to analyze pbmc10k scRNA dataset downloaded from 
`10x official website <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3>`_ (chemistry v3.1). This dataset
has also been used as the demo query data in `Azimuth <https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC>`_. It contains 11,769 single 
cells before filtering.

Here we first conduct basic single cell analysis to obtain Leiden clustering results, however, at various resolutions (r=1,2,3). Smaller resolutions lead to
broader clusters, and larger resolution value will result in more granular clustering. We leverage scTriangulate to take the three resolutions as the query 
annotations, and automatically mix-and-match cluster boundary from different resolutions, which at the end, yield scTriangulate reconciled cluster solutions.

Download and preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First load the packages::

    import os
    import sys
    import scanpy as sc
    from sctriangulate import *
    from sctriangulate.preprocessing import *

The h5 file can be downloaded from `here <https://drive.google.com/file/d/1_s-a621ADH5Y3cHW32WusFyoOo5cKs_b/view?usp=sharing>`_. We used scanpy and scTriangulate
preprocessing module to conduct basic QC filtering and single cell pipeline::

    adata = sc.read_10x_h5('./pbmc_10k_v3_filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

Visualize the important QC metrics and make the decision on the proper cutoffs::

    for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
        sc.pl.violin(adata,key,jitter=0.4)
        plt.savefig('qc_violin_{}.pdf'.format(key),bbox_inches='tight')
        plt.close()

    sc.pl.scatter(adata,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
    plt.savefig('qc_scatter.pdf',bbox_inches='tight')
    plt.close()

.. image:: ./_static/tutorial/single_modality/qc_total_counts.png
   :height: 250px
   :width: 320px
   :align: left
   :target: target

.. image:: ./_static/tutorial/single_modality/qc_n_genes_by_counts.png
   :height: 250px
   :width: 320px
   :align: right
   :target: target

.. image:: ./_static/tutorial/single_modality/qc_pct_counts_mt.png
   :height: 250px
   :width: 320px
   :align: left
   :target: target

.. image:: ./_static/tutorial/single_modality/qc_scatter.png
   :height: 250px
   :width: 320px
   :align: right
   :target: target

We filtered out the cells whose min_genes = 300, min_counts = 500, mt > 20%, 11,022 cells left::

    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.filter_cells(adata, min_counts=500)
    adata = adata[adata.obs.pct_counts_mt < 20, :]  
    print(adata)  # 11022 × 33538


Then we will use scTriangulate wrapper functions to obtain the Leiden clutser results at different resolutions (r=1,2,3), specifically, 
we chose number of PCs as 50, and 3000 highly variable genes::

    adata = scanpy_recipe(adata,is_log=False,resolutions=[1,2,3],pca_n_comps=50,n_top_genes=3000)

After running this command, we will have three columns in ``adata.obs``, namely, ``sctri_rna_leiden_1``, ``sctri_rna_leiden_2``, ``sctri_rna_leiden_3``. 
Also a h5ad file named ``adata_after_scanpy_recipe_rna_1_2_3_umap_True.h5ad`` will be automatically saved to current directory so there's no need to re-run this
step again, Now let's visualize them::

    umap_dual_view_save(adata,cols=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
    # three umaps will be saved to your current directory.

.. image:: ./_static/tutorial/single_modality/three_resolutions.png
   :height: 300px
   :width: 900px
   :align: center
   :target: target

As we can see, different resolutions lead to various number of clusters, and it is clear that certain regions got sub-divided in higher resolutions. However,
we don't know whether this sub-populations are valid off the top of our heads. **Here comes scTriangulate, which will scan each clusters at each resolutions,
and mix-and-match different solutions to achieve an optimal one.**

Running scTriangulate
~~~~~~~~~~~~~~~~~~~~~~~~~

lazy_run
+++++++++++

Running scTriangulate can be as simple as two steps, we first instantiate the ``ScTriangulate`` object, then call ``lazy_run`` class function which will
handle every thing for us::

    adata = sc.read('adata_after_scanpy_recipe_rna_1_2_3_umap_True.h5ad')
    sctri = ScTriangulate(dir='./output',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
    sctri.lazy_run()  # done!!!

However for the purpose of instructing users how to understand this tool, we are going to run it step by step. We first instantiate ``ScTriangulate`` object 
by specify:

1. ``dir``, where all the intermediate and final results/plots will go into?
2. ``adata``, the adata that we want to start with.
3. ``query``, a list contains all the annotations that we want to triangulate.

The ``dir`` doesn't need to be an existing folder, the program will automatically create one if not present.

Step1: compute_metrics
+++++++++++++++++++++++++

The first step of running scTriangulate is to calculate the biologically meaningful metrics for each cluster in each resolution, by default, scTriangulate will
use ``reassign score``, ``TFIDF10 score``, ``TFIDF5 score`` and ``SCCAF score`` to measure the robustness and stability of each cluster, the metrics can be modified
through ``sctri.metrics`` attribute list::

    adata = sc.read('adata_after_scanpy_recipe_rna_1_2_3_umap_True.h5ad')
    sctri = ScTriangulate(dir='./output',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
    sctri.compute_metrics(parallel=True,scale_sccaf=True)
    sctri.serialize('break_point_after_metrics.p')   # save it for next step

After this step, 3 * 4 = 12 columns will be added to the ``sctri.adata.obs`` dataframe, 3 means 3 resolutions, 4 means 4 metrics, those columns store the metrics
we just calculated, click `this <https://docs.google.com/spreadsheets/d/1tuaX09ZaYCWAPa4Nq9HQZDspW2lWaU4-/edit#gid=116886659>`_ to have a look at the resultant data frame.

Step2: compute_shapley
++++++++++++++++++++++++

The second step is to utilize the calculated metrics, and assess which annotation/cluster is the best for **each single cell**. So the program iterate each row,
which is a single cell, retrive all the metrics associated with each cluster, and calculate shapley value of each cluster (in this case, each single cell has 
three conflicting clusters). Then the program will assign the cell to the "best" clusters amongst all solutions. We refer the resultant cluster assignment as
``raw`` cluster result::

    sctri = ScTriangulate.deserialize('output/break_point_after_metrics.p')
    sctri.compute_shapley(parallel=True)
    sctri.serialize('break_point_after_shapley.p')

After this step, 3 + 1 + 1 + 1 columns will be added to the ``sctri.adata.obs``, they are 3 columns corresponding to the shapley value for each annotation, plus
one column named 'final_annotation' storing which annotation is the winner for each cell, and column 'raw' contains raw clusters which are basically annotation
name and cluster name but concatenated by `@` symbol. Last added column is 'prefix', which is just a concatenation of original cluster and current raw cluster. `Click
to have a look at the data frame <https://docs.google.com/spreadsheets/d/19HgWRNdOjn087f_7OTrEf_84IU-fgCGE/edit#gid=431987812>`_.

Step3: prune_result
++++++++++++++++++++++++

This step is to prune the raw result, we first evaluate the robustness of the raw clusters using same set of stability metrics and add the relatively unstable
clusters to ``invalid`` category. (win_fraction < 0.25 by default, meaning if a cluster originally has 100 cells, but has only <25 cells left). The cells in these
unstable invalid clusters will be reassigned to its nearest neightbor's cluster label. After this step, we have ``pruned`` reusult::

    sctri = ScTriangulate.deserialize('output/break_point_after_shapley.p')
    sctri.prune_result()
    sctri.serialize('break_point_after_prune.p')


Step4: building the viewer
++++++++++++++++++++++++++++++

We provide an automatically generated webpage, called scTriangulate viewer, to allow users to dynamically navigate the robustness of each cluster from each
annotations (cluster viewer). Also, it enables the inspection of further heterogeneity that might not have been captured by a 
single annotation (hetergeneity viewer). The logics of following codes are simple, we first build html, then we generate the figures that the html page would 
need to render it::

    sctri = ScTriangulate.deserialize('break_point_after_prune.p')
    sctri.viewer_cluster_feature_html()
    sctri.viewer_cluster_feature_figure(parallel=False,select_keys=['sctri_rna_leiden_1','pruned'])
    sctri.viewer_heterogeneity_html(key='sctri_rna_leiden_1')
    sctri.viewer_heterogeneity_figure(key='sctri_rna_leiden_1')



Inspect the results
~~~~~~~~~~~~~~~~~~~~~~

Now we start to look at the scTriangulate results,

Comparison with Azimuth mapping
++++++++++++++++++++++++++++++++++

Azimuth leverages > 200 ADTs to delineate the major populations in PBMC, which can serve as a silver standard. First we obtain the Azimuth mapping results 
using the h5ad object after we performed qc::

    sctri = ScTriangulate.deserialize('output/break_point_after_prune.p')
    add_azimuth(sctri.adata,'azimuth_pred.tsv')
    for col in ['azimuth','pruned','final_annotation']:
        sctri.plot_umap(col,'category')

.. image:: ./_static/tutorial/single_modality/azimuth.png
   :height: 400px
   :width: 500px
   :align: center
   :target: target

.. image:: ./_static/tutorial/single_modality/final_annotation.png
   :height: 400px
   :width: 500px
   :align: center
   :target: target

.. image:: ./_static/tutorial/single_modality/pruned.png
   :height: 400px
   :width: 500px
   :align: center
   :target: target

As you can see, scTriangulate can mix-and-match different resolutions, shown in the ``final_annotation`` column, and the merged final results have good 
agreement with Azimuth. 

Discover hidden heterogeneity
+++++++++++++++++++++++++++++++++

scTrangulate, by design, could greedily discover any hidden heterogeneity via levaraging the cluster boundaries from each annotation. Here the scTriangulate 
suggests sub-dividing of CD14 Mono population which has been annotated in Azimuth reference.














Multi-modal workflow
-----------------------------------

