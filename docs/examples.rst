Examples
==========

Different resolutions
------------------------

We are going to demonstrate how to reconcile different resolutions when using unsupervised clustering algorithm (i.e. Leiden).

loading dataset::

    sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')


Different references
------------------------

We are going to demonstrate how to reconcile different mapping results from diverse reference atlas. 

loading dataset::

    sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')

Different modalities
-----------------------------------

We are going to demonstrate how to reconcile different clustering results from matched multi-modal measurements.

**1. CITE-Seq (RNA + ADT)**

loading dataset::

    sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')

**2. Multiome (RNA + ATAC)**

loading dataset::

    sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')

**3. TEA-Seq (RNA + ADT + ATAC)**

loading dataset::

    sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')