API
======

ScTriangulate Class Methods
-----------------------------

.. autoclass:: sctriangulate.main_class.ScTriangulate
    :members: 
    :exclude-members: confusion_to_df, plot_heterogeneity, gene_to_df, get_metrics_and_shapley, 
                      salvage_run, lazy_run, add_to_invalid, add_to_invalid_by_win_fraction, clear_invalid,
                      serialize, deserialize, add_new_metrics, plot_winners_statistics, plot_clusterability,
                      display_hierarchy, doublet_predict, compute_metrics, run_single_key_assessment, penalize_artifact,
                      regress_out_size_effect, compute_shapley, pruning, plot_umap, plot_confusion, plot_cluster_feature,
                      modality_contributions, plot_multi_modal_feature_rank, plot_long_heatmap, viewer_cluster_feature_figure,
                      viewer_cluster_feature_html, viewer_heterogeneity_figure, viewer_heterogeneity_html, plot_concordance


(static) salvage_run()
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.salvage_run

(statis) deserialize()
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.deserialize

add_new_metrics()
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.add_new_metrics

add_to_invalid()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.add_to_invalid

add_to_invalid_by_win_fraction()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.add_to_invalid_by_win_fraction

clear_invalid()
~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.clear_invalid

compute_metrics()
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.compute_metrics

compute_shapley()
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.compute_shapley

confusion_to_df()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.confusion_to_df

display_hierarchy()
~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.display_hierarchy

doublet_predict()
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.doublet_predict

gene_to_df()
~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.gene_to_df

get_metrics_and_shapley()
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.get_metrics_and_shapley


lazy_run()
~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.lazy_run

modality_contributions()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.modality_contributions

penalize_artifact()
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.penalize_artifact

plot_clusterability()
~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_clusterability

plot_cluster_feature()
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_cluster_feature

plot_concordance()
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_concordance

plot_confusion()
~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_confusion

plot_heterogeneity()
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_heterogeneity

plot_long_heatmap()
~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_long_heatmap


plot_multi_modal_feature_rank()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_multi_modal_feature_rank

plot_umap()
~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_umap

plot_winners_statistics()
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.plot_winners_statistics

pruning()
~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.pruning

regress_out_size_effect()
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.regress_out_size_effect

run_single_key_assessment()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.run_single_key_assessment

serialize()
~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.serialize


viewer_cluster_feature_figure()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.viewer_cluster_feature_figure

viewer_cluster_feature_html()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.viewer_cluster_feature_html

viewer_heterogeneity_figure()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.viewer_heterogeneity_figure

viewer_heterogeneity_html()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.main_class.ScTriangulate.viewer_heterogeneity_html


Preprocessing Module
----------------------

GeneConvert Class
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: sctriangulate.preprocessing.GeneConvert
    :members: 

Normalization Class
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: sctriangulate.preprocessing.Normalization
    :members: 

small_txt_to_adata()
~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.small_txt_to_adata

large_txt_to_mtx()
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.large_txt_to_mtx

mtx_to_adata()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.mtx_to_adata

mtx_to_large_txt()
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.mtx_to_large_txt


add_azimuth()
~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.add_azimuth

add_annotations()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.add_annotations

.. _reference_to_add_umap:

add_umap()
~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.add_umap

doublet_predict()
~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.doublet_predict

make_sure_adata_writable()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.make_sure_adata_writable

scanpy_recipe()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.scanpy_recipe

concat_rna_and_other()
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.concat_rna_and_other

umap_dual_view_save()
~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.umap_dual_view_save

make_sure_mat_dense()
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.make_sure_mat_dense

make_sure_mat_sparse()
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.make_sure_mat_sparse

gene_activity_count_matrix_old_10x()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.gene_activity_count_matrix_old_10x

gene_activity_count_matrix_new_10x()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.gene_activity_count_matrix_new_10x

find_genes()
~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.find_genes

reformat_peak()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.preprocessing.reformat_peak


Colors Module
----------------------

retrieve_pretty_cmap()
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.retrieve_pretty_cmap

retrieve_pretty_colors()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.retrieve_pretty_colors

pick_n_colors()
~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.pick_n_colors

colors_for_set()
~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.colors_for_set

bg_greyed_cmap()
~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.bg_greyed_cmap

generate_block()
~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.generate_block

generate_gradient()
~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: sctriangulate.colors.generate_gradient



