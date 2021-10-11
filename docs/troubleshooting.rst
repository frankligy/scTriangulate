Troubleshooting
==================

Error getting the Enrichr libraries
-------------------------------------------------

Error like that::

    Traceback (most recent call last):
    File "./run.py", line 18, in <module>
    sctri.lazy_run(scale_sccaf=False)
    File "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/pbmc3k_test/sctriangulate/main_class.py", line 383, in lazy_run
    self.compute_metrics(parallel=compute_metrics_parallel,scale_sccaf=scale_sccaf)
    File "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/pbmc3k_test/sctriangulate/main_class.py", line 704, in compute_metrics
    collect = collect.get()
    File "/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/lib/python3.6/multiprocessing/pool.py", line 644, in get
    raise self._value
    Exception: Error getting the Enrichr libraries

The reason is, the GESAPY package, which performs gene enrichment analysis automatically in the scTriangulate program, needs internet connection.
So if you are using Linux high-performance compute environment, please make sure the program can have access to Enrichr database. This restriction 
maybe resolved in the future version.


Job killed because reaching the maximum RAM
------------------------------------------------

The reason is, by default scTriangulate will utilize python standard library multiprocessing to automatically employ multiple cores for faster
computation, however, it comes at the cost of memory footprint. In the cases of very large dataset (100k cells or 100k features) with more than 5 
annotation-sets, it may reach the maximum RAM that you allocate. If that happen, just simply change the argument ``compute_metrics_parallel=False``::

  sctri.lazy_run(compute_metrics_parallel=False)






