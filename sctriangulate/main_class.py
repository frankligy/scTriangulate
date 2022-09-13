import sys
import os
import copy
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import to_rgba
import matplotlib as mpl
import seaborn as sns
from anytree import Node, RenderTree
from scipy.sparse import issparse,csr_matrix
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage,leaves_list
import multiprocessing as mp
import platform
import subprocess
import re

import scanpy as sc
import anndata as ad
import gseapy as gp
import scrublet as scr

from .logger import *
from .shapley import *
from .metrics import *
from .viewer import *
from .prune import *
from .colors import *
from .preprocessing import *


import matplotlib as mpl

# for publication ready figure
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'



def sctriangulate_setting(backend='Agg',png=False):
    # change the backend
    mpl.use(backend)
    if png:
        # for publication and super large dataset
        mpl.rcParams['savefig.dpi'] = 600
        mpl.rcParams['figure.dpi'] = 600


# define ScTriangulate Object
class ScTriangulate(object):

    '''
    How to create/instantiate ScTriangulate object.

    :param dir: Output folder path on the disk, will create if not exist
    :param adata: input adata file
    :param query: a python list contains the annotation names to query
    :param species: string, either human (default) or mouse, it will impact how the program searches for artifact genes in the database
    :param criterion: int, it controls what genes would be considered as artifact genes:

            * `criterion1`: all will be artifact
            * `criterion2`: all will be artifact except cellcycle  [Default]
            * `criterion3`: all will be artifact except cellcycle, ribosome
            * `criterion4`: all will be artifact except cellcycle, ribosome, mitochondrial
            * `criterion5`: all will be artifact except cellcycle, ribosome, mitochondrial, antisense
            * `criterion6`: all will be artifact except cellcycle, ribosome, mitochondrial, antisense, predict_gene

    :param verbose: int, it controls how the log file will be generated. 1 means print to stdout (default), 2 means print to a file in the directory
                    specified by dir parameter.
    :param add_metrics: python dictionary. These allows users to add additional metrics to favor or disqualify certain cluster. By default,
                        we add tfidf5 score {'tfidf5':tf_idf5_for_cluster}, remember the value in the dictionary should be the name of a callable, user
                        can define the callable by themselves. If don't want any addded metrics, using empty dict {}.
    .. note:

        For the callable, the signature should be func(adata,key,**kwargs) -> mapping {cluster1:0.5,cluster2:0.6}, when running the program in
        lazy_run function, we need to specify added_metrics_kwargs as a list, each element in the list is a dictionary that corresponds to the kwargs
        that will be passed to each callable. 

    :param predict_doublet: boolean or string, whether to predict doublet using scrublet or not. Valid value:

        * True: will predict doublet score
        * False: will not predict doublet score
        * (string) precomputed: will not predict doublet score but just use existing one

    Example::

        adata = sc.read('pbmc3k_azimuth_umap.h5ad')
        sctri = ScTriangulate(dir='./output',adata=adata,query=['leiden1','leiden2','leiden3'])

    '''

    def __init__(self,dir,adata,query,species='human',criterion=2,verbose=1,reference=None,add_metrics={'tfidf5':tf_idf5_for_cluster},
                    predict_doublet=False):

        self.verbose = verbose
        self.dir = dir
        self._create_dir_if_not_exist()
        self.adata = adata
        self.query = query
        if reference is None:
            self.reference = self.query[0]
        else:
            self.reference = reference
        self.species = species
        self.criterion = criterion
        self.score = {}
        self.cluster = {}
        self.uns = {}
        self.metrics = ['reassign','tfidf10','SCCAF','doublet']   # default metrics
        self.add_metrics = {}                         # user can add their own, key is metric name, value is callable
        self.total_metrics = self.metrics               # all metrics considered

        self._set_logging()          
        self._check_adata()
        self.size_dict, _ = get_size(self.adata.obs,self.query)
        self.invalid = []

        # run doublet predict by default in the initialization
        if predict_doublet:
            if not predict_doublet == 'precomputed':
                self.doublet_predict()
        else:
            logger_sctriangulate.info('skip scrublet doublet prediction, instead doublet is filled using value 0.5')
            doublet_scores = np.full(shape=self.adata.obs.shape[0],fill_value=0.5)  # add a dummy score
            self.adata.obs['doublet_scores'] = doublet_scores

        # add add_metrics by default in the initialization
        self.add_new_metrics(add_metrics)


    def __str__(self):  # when you print(instance) in REPL
        return 'ScTriangualate Object:\nWorking directory is {0}\nQuery Annotation: {1}\nReference Annotation: {2}\n'\
            'Species: {3}\nCriterion: {4}\nTotal Metrics: {5}\nScore slot contains: {6}\nCluster slot contains: {7}\nUns slot contains: {8}\n'\
            'Invalid cluster: {9}'.format(self.dir, self.query,self.reference,self.species,self.criterion,self.total_metrics, list(self.score.keys()),
            list(self.cluster.keys()),list(self.uns.keys()),self.invalid)

    def __repr__(self):  # when you type the instance in REPL
        return 'ScTriangualate Object:\nWorking directory is {0}\nQuery Annotation: {1}\nReference Annotation: {2}\n'\
            'Species: {3}\nCriterion: {4}\nTotal Metrics: {5}\nScore slot contains: {6}\nCluster slot contains: {7}\nUns slot contains: {8}\n'\
            'Invalid cluster: {9}'.format(self.dir, self.query,self.reference,self.species,self.criterion,self.total_metrics, list(self.score.keys()),
            list(self.cluster.keys()),list(self.uns.keys()),self.invalid)

    def _create_dir_if_not_exist(self):
        if not os.path.exists(self.dir):
            os.mkdir(self.dir)

    def _check_adata(self):
        # step1: make all cluster name str
        if self.reference in self.query:
            all_keys = self.query
        else:
            all_keys = copy.deepcopy(self.query)
            all_keys.append(self.reference)
        for key in all_keys:
            self.adata.obs[key] = self.adata.obs[key].astype('str')
            self.adata.obs[key] = self.adata.obs[key].astype('category')
        # step2: replace invalid char in cluster and key name    
        ## replace cluster name
        invalid_chars = ['/','@','$',' ']
        for key in all_keys:
            for ichar in invalid_chars:
                self.adata.obs[key] = self.adata.obs[key].str.replace(ichar,'_')
        
        ## replace key name
        for key in all_keys:
            for ichar in invalid_chars:
                self.adata.obs.rename(columns={key:key.replace(ichar,'_')},inplace=True)   

        ## replace query as well
        tmp = []
        for item in self.query:
            for ichar in invalid_chars:
                item = item.replace(ichar,'_')
            tmp.append(item)
        self.query = tmp

        ## replace reference as well
        new = self.reference
        for ichar in invalid_chars:
            new = new.replace(ichar,'_')
        self.reference = new

        # step3: remove index name for smooth h5ad writing
        self.adata.obs.index.name = None
        self.adata.var.index.name = None

    def _set_logging(self):
        import warnings
        warnings.simplefilter("ignore")

        # get all logger and make them silent
        for pkg in ['scanpy','gseapy','scrublet']:
            logging.getLogger(pkg).setLevel(logging.CRITICAL)

        # configure own logger
        if self.verbose == 1:
            c_handler = logging.StreamHandler()
            c_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s' )
            c_handler.setFormatter(c_formatter)        
            logger_sctriangulate.addHandler(c_handler)
            logger_sctriangulate.setLevel(logging.INFO)  # you can not setLevel for the Handler, as the root logger already has default handler and level is Warning, just add handler will not gonna change anything
            logger_sctriangulate.info('Choosing logging to console (VERBOSE=1)')

        elif self.verbose == 2:
            if not os.path.exists(self.dir):
                os.mkdir(self.dir)
            f_handler = logging.FileHandler(os.path.join(self.dir,'scTriangulate.log'))
            f_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s' )
            f_handler.setFormatter(f_formatter)
            logger_sctriangulate.addHandler(f_handler)
            logger_sctriangulate.setLevel(logging.INFO)
            logger_sctriangulate.info('Choosing logging to a log file (VERBOSE=2)')


    def _to_dense(self):
        self.adata.X = self.adata.X.toarray() 
        
    def _to_sparse(self):
        self.adata.X = csr_matrix(self.adata.X)

    def obs_to_df(self,name='sctri_inspect_obs.txt'):
        self.adata.obs.to_csv(os.path.join(self.dir,name),sep='\t')

    def var_to_df(self,name='sctri_inspect_var.txt'):
        self.adata.var.to_csv(os.path.join(self.dir,name),sep='\t')

    def gene_to_df(self,mode,key,raw=False,col='purify',n=100):
        '''
        Output {mode} genes for all clusters in one annotation (key), mode can be either 'marker_genes' or 'exclusive_genes'.

        :param mode: python string, either 'marker_genes' or 'exclusive_genes'
        :param key: python string, annotation name
        :param raw: False will generate non-raw (human readable) format. Default: False
        :param col: Only when mode=='marker_genes', whether output 'whole' column or 'purify' column. Default: purify
        :param n: Only when mode=='exclusive_genes', how many top exclusively expressed genes will be printed for each cluster.

        Examples::

            sctri.gene_to_df(mode='marker_genes',key='annotation1')
            sctri.gene_to_df(mode='exclusive_genes',key='annotation1')
        
        
        '''
        if not raw: # reformat the output to human readable
            df = self.uns['{}'.format(mode)][key]
            if mode == 'marker_genes':
                result = pd.Series()
                for i in range(df.shape[0]):
                    cluster = df.index[i]
                    markers = df.iloc[i][col]
                    single_column = pd.Series(data=markers,name=cluster)
                    result = pd.concat([result,single_column],axis=1,ignore_index=True)
                result.drop(columns=0,inplace=True)
                all_clusters = df.index
                result.columns = all_clusters

            elif mode == 'exclusive_genes':
                result = pd.DataFrame({'cluster':[],'gene':[],'score':[]})
                for i in range(df.shape[0]): #  here the exclusive gene df is actually a series
                    cluster = df.index[i]
                    gene = df[i]
                    col_cluster = np.full(n,fill_value=cluster)
                    col_gene = list(gene.keys())[:n]
                    col_score = list(gene.values())[:n]
                    chunk = pd.DataFrame({'cluster':col_cluster,'gene':col_gene,'score':col_score})
                    result = pd.concat([result,chunk],axis=0)
            result.to_csv(os.path.join(self.dir,'sctri_gene_to_df_{}_{}.txt'.format(mode,key)),sep='\t',index=None)

        elif raw:
            self.uns['{}'.format(mode)][key].to_csv(os.path.join(self.dir,'sctri_gene_to_df_{}_{}.txt'.format(mode,key)),sep='\t')

    def extract_stability(self,keys=None):
        '''
        To extract cluster stability information

        :params keys: a list, containing the annotation column names, None means all in self.query

        Examples::

            sctri.extract_stability(keys=['annotation1','annotation2'])
        '''
        if keys is None:
            keys = self.query
        for key in keys:
            series_list = []
            for metric in self.total_metrics:
                cluster_to_score = self.score[key]['cluster_to_{}'.format(metric)]
                series_list.append(pd.Series(data=cluster_to_score,name=metric))
            pd.concat(series_list,axis=1).to_csv(os.path.join(self.dir,'stability_{}.txt'.format(key)),sep='\t')

    def confusion_to_df(self,mode,key):
        '''
        Print out the confusion matrix with cluster labels (dataframe).

        :param mode: either 'confusion_reassign' or 'confusion_sccaf'
        :param mode: python string, for example, 'annotation1'

        Examples::

            sctri.confusion_to_df(mode='confusion_reassign',key='annotation1')
        
        '''
        self.uns['{}'.format(mode)][key].to_csv(os.path.join(self.dir,'sctri_confusion_to_df_{}_{}.txt'.format(mode,key)),sep='\t')

    def get_metrics_and_shapley(self,barcode,save=True):
        '''
        For one single cell, given barcode/or other unique index, generate the all conflicting cluster from each annotation,
        along with the metrics associated with each cluster, including shapley value.

        :param barcode: string, the barcode for the cell you want to query.
        :param save: save the returned dataframe to directory or not. Default: True
        :return: DataFrame

        Examples::

            sctri.confusion_to_df(barcode='AAACCCACATCCAATG-1',save=True)

        .. image:: ./_static/get_metrics_and_shapley.png
            :height: 100px
            :width: 700px
            :align: center
            :target: target
        '''
        obs = self.adata.obs
        query = self.query
        total_metrics = self.total_metrics
        row = obs.loc[barcode,:]
        metrics_cols = [j + '@' + i for i in query for j in total_metrics]
        shapley_cols = [i + '_' + 'shapley' for i in query]
        row_metrics = row.loc[metrics_cols].values.reshape(len(query),len(total_metrics))
        df = pd.DataFrame(data=row_metrics,index=query,columns=total_metrics)
        row_shapley = row.loc[shapley_cols].values
        df['shapley'] = row_shapley
        row_cluster = row.loc[query].values
        df['cluster'] = row_cluster
        if save:
            df.to_csv(os.path.join(self.dir,'sctri_metrics_and_shapley_df_{}.txt'.format(barcode)),sep='\t')
        return df

    def prune_result(self,win_fraction_cutoff=0.25,reassign_abs_thresh=10,scale_sccaf=False,layer=None,remove1=True,assess_raw=False,added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None}]):
        self.pruning(method='rank',discard=None,scale_sccaf=scale_sccaf,layer=layer,assess_raw=False)
        self.add_to_invalid_by_win_fraction(percent=win_fraction_cutoff)
        self.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=self.reference)
        self.run_single_key_assessment(key='pruned',scale_sccaf=scale_sccaf,layer=layer,added_metrics_kwargs=added_metrics_kwargs)


    @staticmethod
    def salvage_run(step_to_start,last_step_file,outdir=None,scale_sccaf=True,layer=None,added_metrics_kwargs=[{'species':'human','criterion':2,'layer':None}],compute_shapley_parallel=True,
                    shapley_mode='shapley_all_or_none',shapley_bonus=0.01,win_fraction_cutoff=0.25,
                    reassign_abs_thresh=10,assess_raw=False,assess_pruned=True,viewer_cluster=True,viewer_cluster_keys=None,viewer_heterogeneity=True,
                    viewer_heterogeneity_keys=None,nca_embed=False,n_top_genes=3000,other_umap=None,heatmap_scale=None,heatmap_cmap='viridis',heatmap_regex=None,
                    heatmap_direction='include',heatmap_n_genes=None,heatmap_cbar_scale=None):
        '''
        This is a static method, which allows to user to resume running scTriangulate from certain point, instead of running from very 
        beginning if the intermediate files are present and intact.

        :param step_to_start: string, now support 'assess_pruned'.
        :param last_step_file: string, the path to the intermediate from which we start the salvage run.
        :param outdir: None or string, whether to change the outdir or not.
        
        Other parameters are the same as ``lazy_run`` function.

        Examples::

            ScTriangulate.salvage_run(step_to_start='assess_pruned',last_step_file='output/after_rank_pruning.p')

        '''
        # before running this function, make sure previously generated file/folder are renamed, otherwise, they will be overwritten.
        sctri = ScTriangulate.deserialize(last_step_file)
        if outdir is not None:
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            sctri.dir = outdir
            
        if step_to_start == 'assess_pruned':
            
            sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
            sctri.add_to_invalid_by_win_fraction(percent=win_fraction_cutoff)
            sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=sctri.reference)
            sctri.plot_umap('pruned','category')
            if nca_embed:
                adata = nca_embedding(sctri.adata,10,'pruned','umap',n_top_genes=n_top_genes)
                adata.write(os.path.join(sctri.dir,'adata_nca.h5ad'))
            if assess_pruned:
                sctri.run_single_key_assessment(key='pruned',scale_sccaf=scale_sccaf,layer=layer,added_metrics_kwargs=added_metrics_kwargs)
                sctri.serialize(name='after_pruned_assess.p')
            if viewer_cluster:
                sctri.viewer_cluster_feature_html()
                sctri.viewer_cluster_feature_figure(parallel=False,select_keys=viewer_cluster_keys,other_umap=other_umap)
            if viewer_heterogeneity:
                if viewer_heterogeneity_keys is None:
                    viewer_heterogeneity_keys = [sctri.reference]
                for key in viewer_heterogeneity_keys:
                    sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=key)
                    sctri.viewer_heterogeneity_html(key=key)
                    sctri.viewer_heterogeneity_figure(key=key,other_umap=other_umap,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,
                                                      heatmap_direction=heatmap_direction,heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)
        elif step_to_start == 'build_all_viewers':

            if viewer_cluster:
                sctri.viewer_cluster_feature_html()
                sctri.viewer_cluster_feature_figure(parallel=False,select_keys=viewer_cluster_keys,other_umap=other_umap)
            if viewer_heterogeneity:
                if viewer_heterogeneity_keys is None:
                    viewer_heterogeneity_keys = [sctri.reference]
                for key in viewer_heterogeneity_keys:
                    sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=key)
                    sctri.viewer_heterogeneity_html(key=key)
                    sctri.viewer_heterogeneity_figure(key=key,other_umap=other_umap,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,
                                                      heatmap_direction=heatmap_direction,heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)

        elif step_to_start == 'run_shapley':

            sctri.compute_shapley(parallel=compute_shapley_parallel,mode=shapley_mode,bonus=shapley_bonus)
            sctri.serialize(name='after_shapley.p')
            sctri.pruning(method='rank',discard=None,scale_sccaf=scale_sccaf,layer=layer,assess_raw=assess_raw)
            sctri.serialize(name='after_rank_pruning.p')
            sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
            sctri.add_to_invalid_by_win_fraction(percent=win_fraction_cutoff)
            sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=sctri.reference)
            for col in ['final_annotation','pruned']:
                sctri.plot_umap(col,'category')
            if nca_embed:
                logger_sctriangulate.info('starting to do nca embedding')
                adata = nca_embedding(sctri.adata,10,'pruned','umap',n_top_genes=3000)
                adata.write(os.path.join(sctri.dir,'adata_nca.h5ad'))
            if assess_pruned:
                sctri.run_single_key_assessment(key='pruned',scale_sccaf=scale_sccaf,layer=layer,added_metrics_kwargs=added_metrics_kwargs)
                sctri.serialize(name='after_pruned_assess.p')
            if viewer_cluster:
                sctri.viewer_cluster_feature_html()
                sctri.viewer_cluster_feature_figure(parallel=False,select_keys=viewer_cluster_keys,other_umap=other_umap)
            if viewer_heterogeneity:
                if viewer_heterogeneity_keys is None:
                    viewer_heterogeneity_keys = [self.reference]
                for key in viewer_heterogeneity_keys:
                    sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=key)
                    sctri.viewer_heterogeneity_html(key=key)
                    sctri.viewer_heterogeneity_figure(key=key,other_umap=other_umap,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,
                                                    heatmap_direction='include',heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)

        elif step_to_start == 'run_pruning':

            sctri.pruning(method='rank',discard=None,scale_sccaf=scale_sccaf,layer=layer,assess_raw=assess_raw)
            sctri.serialize(name='after_rank_pruning.p')
            sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
            sctri.add_to_invalid_by_win_fraction(percent=win_fraction_cutoff)
            sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=sctri.reference)
            for col in ['final_annotation','pruned']:
                sctri.plot_umap(col,'category')
            if nca_embed:
                logger_sctriangulate.info('starting to do nca embedding')
                adata = nca_embedding(sctri.adata,10,'pruned','umap',n_top_genes=3000)
                adata.write(os.path.join(sctri.dir,'adata_nca.h5ad'))
            if assess_pruned:
                sctri.run_single_key_assessment(key='pruned',scale_sccaf=scale_sccaf,layer=layer,added_metrics_kwargs=added_metrics_kwargs)
                sctri.serialize(name='after_pruned_assess.p')
            if viewer_cluster:
                sctri.viewer_cluster_feature_html()
                sctri.viewer_cluster_feature_figure(parallel=False,select_keys=viewer_cluster_keys,other_umap=other_umap)
            if viewer_heterogeneity:
                if viewer_heterogeneity_keys is None:
                    viewer_heterogeneity_keys = [self.reference]
                for key in viewer_heterogeneity_keys:
                    sctri.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=key)
                    sctri.viewer_heterogeneity_html(key=key)
                    sctri.viewer_heterogeneity_figure(key=key,other_umap=other_umap,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,
                                                    heatmap_direction='include',heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)



            

    def lazy_run(self,compute_metrics_parallel=True,scale_sccaf=False,layer=None,cores=None,added_metrics_kwargs=[{'species':'human','criterion':2,'layer':None}],compute_shapley_parallel=True,
                 shapley_mode='shapley_all_or_none',shapley_bonus=0.01,win_fraction_cutoff=0.25,reassign_abs_thresh=10,
                 assess_raw=False,assess_pruned=False,viewer_cluster=False,viewer_cluster_keys=None,viewer_heterogeneity=False,viewer_heterogeneity_keys=None,
                 nca_embed=False,n_top_genes=3000,other_umap=None,heatmap_scale=None,heatmap_cmap='viridis',heatmap_regex=None,heatmap_direction='include',heatmap_n_genes=None,
                 heatmap_cbar_scale=None):
        '''
        This is the highest level wrapper function for running every step in one goal.

        :param compute_metrics_parallel: boolean, whether to parallelize ``compute_metrics`` step. Default: True
        :param scale_sccaf: boolean, whether to first scale the expression matrix before running sccaf score. Default: True
        :param layer: None or str, the adata layer where the raw count is stored, useful when calculating tfidf score when adata.X has been skewed (no zero value, like totalVI denoised value)
        :param cores: None or int, how many cores you'd like to specify, by default, it is min(n_annotations,n_available_cores) for metrics computing, and n_available_cores for other parallelizable operations
        :param added_metrics_kwargs: list, see the notes in __init__ function, this is to specify additional arguments that will be passed to each added metrics callable.
        :param compute_shapley_parallel: boolean, whether to parallelize ``compute_parallel`` step. Default: True
        :param shapley_mode: string, accepted values:

                     * `shapley_all_or_none`: default, computing shapley, and players only get points when it beats all
                     * `shapley`: computing shapley, but players get points based on explicit ranking, say 3 players, if ranked first, you get 3, if running up, you get 2
                     * `rank_all_or_none`: no shapley computing, importance based on ranking, and players only get points when it beats all
                     * `rank`: no shapley computing, importance based on ranking, but players get points based on explicit ranking as described above

        :param shapley_bonus: float, default is 0.01, an offset value so that if the runner up is just {bonus} inferior to first place, it will still be a valid cluster
        :param win_fraction_cutoff: float, between 0-1, the cutoff for function ``add_invalid_by_win_fraction``. Default: 0.25
        :param reassign_abs_thresh: int, the cutoff for minimum number of cells a valid cluster should haves. Default: 10
        :param assess_raw: boolean, whether to run the same cluster assessment metrics on raw cluster labels. Default: False
        :param assess_pruned: boolean, whether to run same cluster assessment metrics on final pruned cluster labels. Default: True
        :param viewer_cluster: boolean, whether to build viewer html page for all clusters' diagnostic information. Default: True
        :param viewer_cluster_keys: list, clusters from what annotations we want to view on the viewer, only clusters within this annotation whose diagnostic
                                    plot will be generated under the dir name *figure4viewer*. Default: None, means all annotations in the sctri.query will be 
                                    used.
        :param viewer_heterogeneity: boolean, whether to build the viewer to show the heterogeneity based on one reference annotation. Default: True
        :param viewer_heterogeneity_keys: list, the annotations we want to serve as the reference. Default: None, means the first annotation in sctri.query
                                          will be used as the reference.
        
        Examples::

            sctri.lazy_run(viewer_heterogeneity_keys=['annotation1','annotation2'])
        '''
        logger_sctriangulate.info('Starting to compute stability metrics, ignore scanpy logging like "Trying to set..." or "Storing ... as categorical"')
        self.compute_metrics(parallel=compute_metrics_parallel,scale_sccaf=scale_sccaf,layer=layer,added_metrics_kwargs=added_metrics_kwargs,cores=cores)
        self.serialize(name='after_metrics.p')
        logger_sctriangulate.info('Starting to compute shapley')
        self.compute_shapley(parallel=compute_shapley_parallel,mode=shapley_mode,bonus=shapley_bonus,cores=cores)
        self.serialize(name='after_shapley.p')
        logger_sctriangulate.info('Starting to prune and reassign the raw result to get pruned results')
        self.pruning(method='rank',discard=None,scale_sccaf=scale_sccaf,layer=layer,assess_raw=assess_raw)
        self.serialize(name='after_rank_pruning.p')
        self.uns['raw_cluster_goodness'].to_csv(os.path.join(self.dir,'raw_cluster_goodness.txt'),sep='\t')
        self.add_to_invalid_by_win_fraction(percent=win_fraction_cutoff)
        self.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=self.reference)
        for col in ['final_annotation','pruned']:
            self.plot_umap(col,'category')
        # output necessary result
        make_sure_adata_writable(self.adata)
        self.adata.write(os.path.join(self.dir,'sctriangulate.h5ad'))
        self.adata.obs.to_csv(os.path.join(self.dir,'sctri_barcode2cellmetadata.txt'),sep='\t')
        pd.DataFrame(data=self.adata.obsm['X_umap'],index=self.adata.obs_names,columns=['umap_x','umap_y']).to_csv(os.path.join(self.dir,'sctri_umap_coord.txt'))
        self.extract_stability(keys=self.query)
        for q in self.query:
            self.gene_to_df(mode='marker_genes',key=q)
            self.gene_to_df(mode='exclusive_genes',key=q)

        if nca_embed:
            logger_sctriangulate.info('starting to do nca embedding')
            adata = nca_embedding(self.adata,10,'pruned','umap',n_top_genes=3000)
            adata.write(os.path.join(self.dir,'adata_nca.h5ad'))
        if assess_pruned:
            logger_sctriangulate.info('starting to get stability metrics on pruned final results')
            self.run_single_key_assessment(key='pruned',scale_sccaf=scale_sccaf,layer=layer,added_metrics_kwargs=added_metrics_kwargs)
            self.serialize(name='after_pruned_assess.p')
            subprocess.run(['rm','-r','{}'.format(os.path.join(self.dir,'scTriangulate_local_mode_enrichr/'))])
            # update the old output 
            make_sure_adata_writable(self.adata)
            self.adata.write(os.path.join(self.dir,'sctriangulate.h5ad'))
            self.adata.obs.to_csv(os.path.join(self.dir,'sctri_barcode2cellmetadata.txt'),sep='\t')
            self.extract_stability(keys=['pruned'])
            self.gene_to_df(mode='marker_genes',key='pruned')
            self.gene_to_df(mode='exclusive_genes',key='pruned')
        if viewer_cluster:
            self.viewer_cluster_feature_html()
            self.viewer_cluster_feature_figure(parallel=False,select_keys=viewer_cluster_keys,other_umap=other_umap)
        if viewer_heterogeneity:
            if viewer_heterogeneity_keys is None:
                viewer_heterogeneity_keys = [self.reference]
            for key in viewer_heterogeneity_keys:
                self.pruning(method='reassign',abs_thresh=reassign_abs_thresh,remove1=True,reference=key)
                self.viewer_heterogeneity_html(key=key)
                self.viewer_heterogeneity_figure(key=key,other_umap=other_umap,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,
                                                 heatmap_direction='include',heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)


            

    def add_to_invalid(self,invalid):
        '''
        add individual raw cluster names to the sctri.invalid attribute list.

        :param invalid: list or string, contains the raw cluster names to add

        Examples::

            sctri.add_to_invalid(invalid=['annotation1@c3','annotation2@4'])
            sctri.add_to_invalid(invalid='annotation1@3')

        '''
        try:
            self.invalid.extend(invalid)
        except AttributeError:
            self.invalid = []
            self.invalid.extend(invalid)
        finally:
            tmp = list(set(self.invalid))
            self.invalid = tmp

    def add_to_invalid_by_win_fraction(self,percent=0.25):
        '''
        add individual raw cluster names to the sctri.invalid attribute list by win_fraction

        :param percent: float, from 0-1, the fraction of cells within a cluster that were kept after the game. Default: 0.25

        Examples::

            sctri.add_to_invalid_by_win_fraction(percent=0.25)
        '''
        df = self.uns['raw_cluster_goodness']
        invalid = df.loc[df['win_fraction']<percent,:].index.tolist()
        self.add_to_invalid(invalid)

    def clear_invalid(self):
        '''
        reset/clear the sctri.invalid to an empty list

        Examples::

            sctri.clear_invalid()
        '''
        del self.invalid
        self.invaild = []

    def serialize(self,name='sctri_pickle.p'):
        '''
        serialize the sctri object through pickle protocol to the disk

        :param name: string, the name of the pickle file on the disk. Default: sctri_pickle.p

        Examples::

            sctri.serialize()
        '''

        with open(os.path.join(self.dir,name),'wb') as f:
            pickle.dump(self,f)

    @staticmethod
    def deserialize(name):
        '''
        This is static method, to deserialize a pickle file on the disk back to the ram as a sctri object

        :param name: string, the name of the pickle file on the disk. 

        Examples::

            ScTriangulate.deserialize(name='after_rank_pruning.p')
        '''
        with open(name,'rb') as f:
            sctri = pickle.load(f)
        sctri._set_logging()
        logger_sctriangulate.info('unpickled {} to memory'.format(name))
        return sctri

    def add_new_metrics(self,add_metrics):
        '''
        Users can add new callable or pre-implemented function to the sctri.metrics attribute.

        :param add_metrics: dictionary like {'metric_name': callable}, the callable can be a string of a scTriangulate pre-implemented function, for example,
                            'tfidf5','tfidf1'. Or a callable.
        
        Examples::

            sctri.add_new_metrics(add_metrics={'tfidf1':tfidf1})  # make sure first from sctriangualte.metrics import tfidf1
        '''
        for metric,func in add_metrics.items():
            self.add_metrics[metric] = func
        self.total_metrics.extend(list(self.add_metrics.keys()))

    def plot_winners_statistics(self,col,fontsize=3,plot=True,save=True):
        '''
        For triangulated clusters, either 'raw' or 'pruned', visualize what fraction of cells won the game.
        A horizontal barplot will be generated and a dataframe with winners statistics will be returned.

        :param col: string, either 'raw' or 'pruned'
        :param fontsize: int, the fontsize for the y-label. Default: 3
        :param plot: boolean, whether to plot or not. Default: True
        :param save: boolean, whether to save the plot to the sctri.dir or not. Default: True

        :return: DataFarme 

        Examples::

            sctri.plot_winners_statistics(col='raw',fontsize=4)

        .. image:: ./_static/plot_winners_statistics.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target

        '''
        new_size_dict = {}  # {gs@ERP4: 100}
        for key,value in self.size_dict.items():
            for sub_key,sub_value in value.items():
                composite_key = key + '@' + sub_key
                composite_value = sub_value
                new_size_dict[composite_key] = composite_value
        obs = self.adata.obs
        winners = obs[col]
        winners_vc = winners.value_counts()
        winners_size = winners_vc.index.to_series().map(new_size_dict).astype('int64')
        winners_prop = winners_vc / winners_size
        winners_stats = pd.concat([winners_vc,winners_size,winners_prop],axis=1)
        winners_stats.columns = ['counts','size','proportion']
        winners_stats.sort_values(by='proportion',inplace=True)
        if plot:
            a = winners_stats['proportion']
            fig,ax = plt.subplots()
            ax.barh(y=np.arange(len(a)),width=[item for item in a.values],color='#FF9A91')
            ax.set_yticks(np.arange(len(a)))
            ax.set_yticklabels([item for item in a.index],fontsize=fontsize)
            ax.set_title('Winners statistics')
            ax.set_xlabel('proportion of cells in each cluster that win')
            if save:
                plt.savefig(os.path.join(self.dir,'winners_statistics.pdf'),bbox_inches='tight')
                plt.close()
        return winners_stats

    def plot_clusterability(self,key,col,fontsize=3,plot=True,save=True):
        '''
        We define clusterability as the number of sub-clusters the program finds out. If a cluster has being suggested to be divided into
        three smaller clusters, then the clueterability of this cluster will be 3.

        :param key: string. The clusters from which annotation that you want to assess clusterability.
        :param col: string. Either 'raw' cluster or 'pruned' cluster.
        :param fontsize: int. The fontsize of x-ticklabels. Default: 3
        :param plot: boolean. Whether to plot the scatterplot or not. Default : True.
        :param save: boolean. Whether to save the plot or not. Default: True

        :return: python dictionary. {cluster1:#sub-clusters}

        Examples::

            sctri.plot_clusterability(key='sctri_rna_leiden_1',col='raw',fontsize=8)

        .. image:: ./_static/plot_clusterability.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target
        

        '''
        bucket = {}   # {ERP4:5}
        obs = self.adata.obs
        for ref,grouped_df in obs.groupby(by=key):
            unique = grouped_df[col].unique()
            bucket[ref] = len(unique)
        bucket = {k: v for k, v in sorted(bucket.items(), key=lambda x: x[1])}
        if plot:
            fig,ax = plt.subplots()
            ax.scatter(x=np.arange(len(bucket)),y=list(bucket.values()),c=pick_n_colors(len(bucket)),s=100)
            ax.set_xticks(np.arange(len(bucket)))
            ax.set_xticklabels(list(bucket.keys()),fontsize=fontsize,rotation=90)
            ax.set_title('{} clusterability'.format(self.reference))
            ax.set_ylabel('clusterabiility: # sub-clusters')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.grid(color='grey',alpha=0.2)
            for i in range(len(bucket)):
                ax.text(x=i,y=list(bucket.values())[i]+0.3,s=list(bucket.keys())[i],ha='center',va='bottom')

            if save:
                plt.savefig(os.path.join(self.dir,'{}_clusterability.pdf'.format(self.reference)),bbox_inches='tight')
                plt.close()
        return bucket


    def display_hierarchy(self,ref_col,query_col,save=True):
        '''
        Display the hierarchy of suggestive sub-clusterings, see the example results down the page.

        :param ref_col: string, the annotation/column name in adata.obs which we want to inspect how it can be sub-divided
        :param query_col: string, any cluster annotation column name
        :param save: boolean, whether to save it to a file or stdout. Default: True

        Examples::

            sctri.display_hierarchy(ref_col='sctri_rna_leiden_1',query_col='raw')

        .. image:: ./_static/display_hierarchy.png
            :height: 400px
            :width: 300px
            :align: center
            :target: target        

        '''
        obs = self.adata.obs
        root = Node(ref_col)
        hold_ref_var = {}
        for ref,grouped_df in obs.groupby(by=ref_col):
            ref_display = '{}[#:{}]'.format(ref,grouped_df.shape[0])
            hold_ref_var[ref] = Node(ref_display,parent=root)
            vf = grouped_df[query_col].value_counts()
            unique = vf.index.tolist()
            if len(unique) == 1: # no sub-clusters
                continue
            else:
                hold_cluster_var = {}
                for item in unique:
                    if vf[item] > 0:
                        item_display = '{}[#:{};prop:{}]'.format(item,vf[item],round(vf[item]/grouped_df.shape[0],2))
                        hold_cluster_var[item] = Node(item_display,parent=hold_ref_var[ref])
        if save:
            with open(os.path.join(self.dir,'display_hierarchy_{}_{}.txt'.format(ref_col,query_col)),'a') as f:
                for pre, fill, node in RenderTree(root):
                    print("%s%s" % (pre, node.name),file=f)
        else:
            for pre, fill, node in RenderTree(root):
                print("%s%s" % (pre, node.name))


    def doublet_predict(self):
        '''
        wrapper function of running scrublet, will add a column on adata.obs called 'doublet_scores'

        Examples::

            sctri.doublet_predict()

        '''
        logger_sctriangulate.info('Running scrublet to get doublet scores, it will take a while and please follow their prompts below:')
        if not issparse(self.adata.X):
            self._to_sparse()
        counts_matrix = self.adata.X  # I don't want adata.X to be modified, so make a copy 
        scrub = scr.Scrublet(counts_matrix)
        doublet_scores,predicted_doublets = scrub.scrub_doublets(min_counts=1,min_cells=1)
        self.adata.obs['doublet_scores'] = doublet_scores
        del counts_matrix
        del scrub
  

    def _add_to_uns(self,name,key,collect):
        try:
            self.uns[name][key] = collect[name]
        except KeyError:
            self.uns[name] = {}
            self.uns[name][key] = collect[name]

    def cluster_performance(self,cluster,competitors,reference,show_cluster_number=False,metrics=None,ylim=None,save=True,format='pdf'):
        '''
        automatic benchmark of scTriangulate clusters with all the individual or competitor annotation, against a 'gold standard' annotation,
        measured by all unsupervised cluster metrics (homogeneity, completeness, v_measure, ARI, NMI).

        :param cluster: string, the scTriangulate annotation column name, for example, pruned
        :param competitiors: list of string, each is a column name of a competitor annotation
        :param reference: string, the column name containing reference annotation, for example, azimuth
        :param show_cluster_number: bool, whether to show the number of cluster of each annotation in the performance line plot
        :param metrics: None or any other, if not None, ARI and NMI will also be plotted
        :param ylim: None or a tuple, specifiying the ylims of plot
        :param save: bool, whether to save the figure
        :param format: string, default is pdf, the format to save

        Examples::

            sctri.cluster_performance(cluster='pruned',competitors=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],
                                      reference='azimuth',show_cluster_number=True,metrics=None)


        .. image:: ./_static/cluster_performance.png
            :height: 350px
            :width: 600px
            :align: center
            :target: target          

        '''
        from sklearn.preprocessing import LabelEncoder
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure, adjusted_mutual_info_score
        result = self.adata.obs
        # label encoder
        reference_encoded = LabelEncoder().fit_transform(result[reference].values)
        competitors_encoded = [LabelEncoder().fit_transform(result[anno].values) for anno in competitors]
        cluster_encoded = LabelEncoder().fit_transform(result[cluster].values)
        # compute metrics for competitors
        ari = []
        ami = []
        homogeneity = []
        completeness = []
        vmeasure = []
        for anno_encoded in competitors_encoded:
            ari.append(adjusted_rand_score(reference_encoded,anno_encoded))
            ami.append(adjusted_mutual_info_score(reference_encoded,anno_encoded))
            h,c,v = homogeneity_completeness_v_measure(reference_encoded,anno_encoded)
            homogeneity.append(h)
            completeness.append(c)
            vmeasure.append(v)
        # compute metrics for cluster
        ari.append(adjusted_rand_score(reference_encoded,cluster_encoded))
        ami.append(adjusted_mutual_info_score(reference_encoded,cluster_encoded))
        h,c,v = homogeneity_completeness_v_measure(reference_encoded,cluster_encoded)
        homogeneity.append(h)
        completeness.append(c)
        vmeasure.append(v)
        # now plot
        fig,ax = plt.subplots()
        ax.plot(np.arange(len(competitors)+1),homogeneity,label='Homogeneity',marker='o',linestyle='--')
        ax.plot(np.arange(len(competitors)+1),completeness,label='Completeness',marker='o',linestyle='--')
        ax.plot(np.arange(len(competitors)+1),vmeasure,label='VMeasure',marker='o',linestyle='--')
        if metrics is not None:
            ax.plot(np.arange(len(competitors)+1),ari,label='ARI',marker='o',linestyle='--')
            ax.plot(np.arange(len(competitors)+1),ami,label='AMI',marker='o',linestyle='--')
        ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
        ax.set_xticks(np.arange(len(competitors)+1))
        ax.set_xticklabels(competitors+[cluster],fontsize=3)
        ax.set_ylabel('Agreement with {}'.format(reference))
        if ylim is not None:
            ax.set_ylim(ylim)
        if show_cluster_number:  # show how many clusters in each annotation
            number = []
            for anno in competitors + [cluster]:
                number.append(len(result[anno].value_counts()))
            for i,num in enumerate(number):
                ax.text(x=i,y=vmeasure[i]+0.01,s=num,fontsize=6)

        if save:
            plt.savefig(os.path.join(self.dir,'cluster_performance_plot.{}'.format(format)),bbox_inches='tight')
            plt.close()

        # assemble returned metrics
        df = pd.DataFrame.from_records(data=[ari,ami,homogeneity,completeness,vmeasure],columns=competitors+[cluster],index=['ari','ami','homogeneity','completeness','vmeasure'])
        if save:
            df.to_csv(os.path.join(self.dir,'cluster_performance.txt'),sep='\t')
        return df


    def compute_metrics(self,parallel=True,scale_sccaf=False,layer=None,added_metrics_kwargs=[{'species': 'human', 'criterion': 2, 'layer': None}],cores=None):
        '''
        main function for computing the metrics (defined by self.metrics) of each clusters in each annotation.
        After the run, q (# query) * m (# metrics) columns will be added to the adata.obs, the column like will be like
        {metric_name}@{query_annotation_name}, i.e. reassign@sctri_rna_leiden_1

        :param parallel: boolean, whether to run in parallel. Since computing metrics for each query annotation is idependent,
                         the program will automatically employ q (# query) cores under the hood. If you want to fully leverage
                         this feature, please make sure you specify at least q (# query) cores when running the program. It is highly
                         recommend to run this in parallel. However, when the dataset is super large and have > 10 query annotation, 
                         we may encounter RAM overhead, in this case, sequential mode will be needed. Default: True
        :param scale_sccaf: boolean, when running SCCAF score, since it is a logistic regression problem at its core, this parameter 
                            controls whether scale the expression data or not. It is recommended to scale the data for any machine learning
                            algorithm, however, the liblinaer solver has been demonstrated to be robust to the scale/unscale options. When
                            the dataset > 50,000 cells or features > 1000000 (ATAC peaks), it is advised to not scale it for faster running time.
        :param layer: None or str, the adata layer where the raw count is stored, useful when calculating tfidf score when adata.X has been skewed (no zero value, 
                      like totalVI denoised value)
        :param added_metrics_kwargs: list, see the notes in __init__ function, this is to specify additional arguments that will be passed to each added metrics callable.
        :param cores: None or int, how many cores you’d like to specify, by default, it is min(n_annotations,n_available_cores) for metrics computing, 
                      and n_available_cores for other parallelizable operations
        
        Examples::

            sctri.compute_metrics(parallel=False)


        '''
        if parallel:
            cores1 = len(self.query)  # make sure to request same numeber of cores as the length of query list
            cores2 = mp.cpu_count()
            cores3 = cores
            if cores is not None:
                cores = cores3
            else:
                cores = min(cores1,cores2)
            logger_sctriangulate.info('Spawn to {} processes'.format(cores))
            pool = mp.Pool(processes=cores)
            self._to_sparse()
            raw_results = [pool.apply_async(each_key_run,args=(self,key,scale_sccaf,layer,added_metrics_kwargs)) for key in self.query]
            pool.close()
            pool.join()
            for collect in raw_results:
                collect = collect.get()
                key = collect['key']
                for metric in self.total_metrics:
                    self.adata.obs['{}@{}'.format(metric,key)] = collect['col_{}'.format(metric)]
                self.score[key] = collect['score_info']
                self.cluster[key] = collect['cluster_info']  
                self._add_to_uns('confusion_reassign',key,collect)
                self._add_to_uns('confusion_sccaf',key,collect)
                self._add_to_uns('marker_genes',key,collect)
                self._add_to_uns('exclusive_genes',key,collect)
            subprocess.run(['rm','-r','{}'.format(os.path.join(self.dir,'scTriangulate_local_mode_enrichr/'))])
            self._to_sparse()

        else:
            logger_sctriangulate.info('choosing to compute metrics sequentially')
            for key in self.query:
                collect = each_key_run(self,key,scale_sccaf,layer,added_metrics_kwargs)
                key = collect['key']
                for metric in self.metrics + list(self.add_metrics.keys()):
                    self.adata.obs['{}@{}'.format(metric,key)] = collect['col_{}'.format(metric)]
                self.score[key] = collect['score_info']
                self.cluster[key] = collect['cluster_info']  
                self._add_to_uns('confusion_reassign',key,collect)
                self._add_to_uns('confusion_sccaf',key,collect)
                self._add_to_uns('marker_genes',key,collect)
                self._add_to_uns('exclusive_genes',key,collect)
            subprocess.run(['rm','-r','{}'.format(os.path.join(self.dir,'scTriangulate_local_mode_enrichr/'))])
            self._to_sparse()

    def run_single_key_assessment(self,key,scale_sccaf,layer,added_metrics_kwargs):
        '''
        this is a very handy function, given a set of annotation, this function allows you to assess the biogical robustness
        based on the metrics we define. The obtained score and cluster information will be automatically saved to self.cluster
        and self.score, and can be further rendered by the scTriangulate viewer.

        :param key: string, the annotation/column name to assess the robustness.
        :scale_sccaf: boolean, whether to scale the expression data before running SCCAF score. See ``compute_metrics`` function
                      for full information.

        Examples::

            sctri.run_single_key_assessment(key='azimuth',scale_sccaf=True)

        '''
        collect = each_key_run(self,key,scale_sccaf,layer,added_metrics_kwargs)
        self._to_sparse()
        self.process_collect_object(collect)

    def process_collect_object(self,collect):
        key = collect['key']
        for metric in self.total_metrics:
            self.adata.obs['{}@{}'.format(metric,key)] = collect['col_{}'.format(metric)]
            self.score[key] = collect['score_info']
            self.cluster[key] = collect['cluster_info']  
            self._add_to_uns('confusion_reassign',key,collect)
            self._add_to_uns('confusion_sccaf',key,collect)
            self._add_to_uns('marker_genes',key,collect)
            self._add_to_uns('exclusive_genes',key,collect)

            

    def penalize_artifact(self,mode,stamps=None,parallel=True):
        '''
        An optional step after running ``compute_metrics`` step and before the ``compute_shapley`` step.
        Basically, we penalize clusters with certain properties by set all their metrics to zero, which forbid them to
        win in the following "game" step. These undesirable properties can be versatial, for example, cellcylce gene 
        enrichment. We current support two mode:

        1. mode1: ``void``, users specifiy which cluster they want to penalize via ``stamps`` parameter.
        2. mode2: ``cellcycle``, program automatically label clusters whose gsea_hit > 5 and gsea_score > 0.8 as invalid cellcyle
           enriched clusters. And those clusters will be penalized.

        :param mode: string, either 'void' or 'cellcycle'.
        :param stamps: list, contains cluster names that the users want to penalize.
        :param parallel: boolean, whether to run this in parallel (scatter and gather). Default: True.

        Examples::

            sctri.penalize_artifact(mode='void',stamps=['sctri_rna_leiden_1@c3','sctri_rna_leiden_2@c5'])
            sctri.penalize_artifact(mode='cellcyle')

        '''

        '''void mode is to set stamp position to 0, stamp is like {leiden1:5}'''
        if mode == 'void':
            obs = self.adata.obs
            self.add_to_invalid(stamps)
            if parallel:
                obs_index = np.arange(obs.shape[0])  # [0,1,2,.....]
                cores = mp.cpu_count()
                sub_indices = np.array_split(obs_index,cores)  # indices for each chunk [(0,1,2...),(56,57,58...),(),....]
                sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_df,sub_df,...]
                pool = mp.Pool(processes=cores)
                logger_sctriangulate.info('spawn {} sub processes for penalizing artifact with mode-{}'.format(cores,mode))
                r = [pool.apply_async(func=penalize_artifact_void,args=(chunk,self.query,stamps,self.total_metrics,)) for chunk in sub_obs]
                pool.close()
                pool.join()
                results = []
                for collect in r:
                    result = collect.get()  # [sub_obs,sub_obs...]
                    results.append(result)
                obs = pd.concat(results)
                self.adata.obs = obs
            else:
                result = penalize_artifact_void(obs,self.query,stamps,self.total_metrics)
                self.adata.obs = result

        elif mode == 'cellcycle':
            # all the clusters that have cell-cycle enrichment > 0 will be collected into stamps
            marker_genes = self.uns['marker_genes']
            stamps = []
            for key,clusters in self.cluster.items():
                for cluster in clusters:
                    gsea_score = marker_genes[key].loc[cluster,:]['gsea']['cellcycle'][0]
                    gsea_hits = marker_genes[key].loc[cluster,:]['gsea']['cellcycle'][1]
                    if gsea_hits > 5 and gsea_score > 0.8:
                        stamps.append(key+'@'+cluster)
            logger_sctriangulate.info('stamps are: {}'.format(str(stamps)))
            self.invalid.extend(stamps)
            obs = self.adata.obs
            if parallel:
                obs_index = np.arange(obs.shape[0])  # [0,1,2,.....]
                cores = mp.cpu_count()
                sub_indices = np.array_split(obs_index,cores)  # indices for each chunk [(0,1,2...),(56,57,58...),(),....]
                sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_df,sub_df,...]
                pool = mp.Pool(processes=cores)
                logger_sctriangulate.info('spawn {} sub processes for penalizing artifact with mode-{}'.format(cores,mode))
                r = [pool.apply_async(func=penalize_artifact_void,args=(chunk,self.query,stamps,self.total_metrics,)) for chunk in sub_obs]
                pool.close()
                pool.join()
                results = []
                for collect in r:
                    result = collect.get()  # [sub_obs,sub_obs...]
                    results.append(result)
                obs = pd.concat(results)
                self.adata.obs = obs
            else:
                result = penalize_artifact_void(obs,self.query,stamps,self.total_metrics)
                self.adata.obs = result
            

    def regress_out_size_effect(self,regressor='background_zscore'):
        '''
        An optional step to regress out potential confounding effect of cluster_size on the metrics. Run after ``compute_metrics`` step
        but before ``compute_shapley`` step. All the metrics in selfadata.obs and self.score will be modified in place.

        :param regressor: string. which regressor to choose, valid values: 'background_zscore', 'background_mean', 'GLM', 'Huber', 'RANSAC', 'TheilSen'

        Example::

            sctri.regress_out_size(regressor='Huber')

        '''

        sctri = self
        '''
        the logic of this function is:
        1, take the score slot of sctriangulate object, reformat to {score:[df_a1,df_a2...],},each df_a is index(c_name),metric,size
        2. for each score, concated df will be subjected to regress_size main function, replace metric in place, deal with NA as well
        3. restore to original score slot {annotation:{score1:{value_dict}}}
        4. map back to each metric column in adata.obs
        '''
        result = {}
        order_of_keys = list(sctri.score.keys())
        for key in sctri.score.keys():
            size = get_size_in_metrics(sctri.adata.obs,key)
            slot = sctri.score[key]
            for score in slot.keys():
                df = pd.concat([pd.Series(slot[score]),pd.Series(size)],axis=1)
                try:
                    result[score].append(df)
                except KeyError:
                    result[score] = []
                    result[score].append(df)

        restore_score = {}
        for key,value in result.items():
            df_inspect_have_na = pd.concat(value,axis=0)
            df_inspect_have_na['ori'] = np.arange(df_inspect_have_na.shape[0])
            mask = df_inspect_have_na[0].isna()
            df_inspect = df_inspect_have_na.dropna(axis=0) # metric, size, ori, index is the cluster names
            df_na = df_inspect_have_na.loc[mask,:] # metric, size, ori, index is the cluster names
            df_inspect[0] = regress_size(df_inspect,regressor=regressor).values # change the metric col to regressed one
            df_na[0] = df_inspect[0].values.min() - 1 # make sure the na has smaller value than non-na ones
            df_all = pd.concat([df_inspect,df_na]).sort_values(by='ori')  # concat and reorder to the original order
            # now need to split it up, back to each annotation df
            rowptr = 0
            chunk_length = [item.shape[0] for item in value]
            for chunkptr,length in enumerate(chunk_length):
                bound = (rowptr,rowptr+length)
                target_df = df_all.iloc[bound[0]:bound[1],:]
                annotation = order_of_keys[chunkptr]
                target_dict = target_df[0].to_dict()
                try:
                    restore_score[annotation][key] = target_dict
                except KeyError:
                    restore_score[annotation] = {}
                    restore_score[annotation][key] = target_dict
                rowptr = bound[1]
        sctri.score = restore_score
        # map all back
        for key in sctri.score.keys():
            for metric in sctri.total_metrics:
                sctri.adata.obs['{}@{}'.format(metric,key)] = sctri.adata.obs[key].map(sctri.score[key]['cluster_to_{}'.format(metric)]).fillna(0).values
        
        return df_inspect_have_na,df_all



    def compute_shapley(self,parallel=True,mode='shapley_all_or_none',bonus=0.01,cores=None):
        '''
        Main core function, after obtaining the metrics for each cluster. For each single cell, let's calculate the shapley
        value for each annotation and assign the cluster to the one with highest shapley value.

        :param parallel: boolean. Whether to run it in parallel. (scatter and gather). Default: True
        :param mode: string, accepted values:

                     * `shapley_all_or_none`: default, computing shapley, and players only get points when it beats all
                     * `shapley`: computing shapley, but players get points based on explicit ranking, say 3 players, if ranked first, you get 3, if running up, you get 2
                     * `rank_all_or_none`: no shapley computing, importance based on ranking, and players only get points when it beats all
                     * `rank`: no shapley computing, importance based on ranking, but players get points based on explicit ranking as described above
        :param bonus: float, default is 0.01, an offset value so that if the runner up is just {bonus} inferior to first place, it will still be a valid cluster
        :param cores: None or int, None will run mp.cpu_counts() to get all available cpus.

        Examples::

            sctri.compute_shapley(parallel=True)

        '''
        if parallel:
            # compute shaley value
            score_colname = copy.deepcopy(self.total_metrics)
            score_colname.remove('doublet')
            data = np.empty([len(self.query),self.adata.obs.shape[0],len(score_colname)])  # store the metric data for each cell
            '''
            data:
            depth is how many sets of annotations
            height is how many cells
            width is how many score metrics
            '''
            for i,key in enumerate(self.query):
                practical_colname = [name + '@' + key for name in score_colname]
                data[i,:,:] = self.adata.obs[practical_colname].values
            final = []
            intermediate = []
            if cores is not None:
                cores = cores
            else:
                cores = mp.cpu_count()
            # split the obs and data, based on cell axis
            obs = self.adata.obs
            obs_index = np.arange(obs.shape[0])
            sub_indices = np.array_split(obs_index,cores)
            sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_obs, sub_obs, sub_obs]
            sub_datas = [data[:,sub_index,:] for sub_index in sub_indices]  # [sub_data,sub_data,....]
            pool = mp.Pool(processes=cores)
            logger_sctriangulate.info('spawn {} sub processes for shapley computing'.format(cores))
            raw_results = [pool.apply_async(func=run_shapley,args=(sub_obs[i],self.query,self.reference,self.size_dict,sub_datas[i],mode,bonus,)) for i in range(len(sub_obs))]
            pool.close()
            pool.join()
            for collect in raw_results: # [(final,intermediate), (), ()...]
                collect = collect.get()
                final.extend(collect[0])
                intermediate.extend(collect[1])
            self.adata.obs['final_annotation'] = final
            decisions = list(zip(*intermediate))
            for i,d in enumerate(decisions):
                self.adata.obs['{}_shapley'.format(self.query[i])] = d

            # get raw sctriangulate result
            obs = self.adata.obs
            obs_index = np.arange(obs.shape[0])  # [0,1,2,.....]
            sub_indices = np.array_split(obs_index,cores)  # indices for each chunk [(0,1,2...),(56,57,58...),(),....]
            sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_df,sub_df,...]
            pool = mp.Pool(processes=cores)
            r = pool.map_async(run_assign,sub_obs)
            pool.close()
            pool.join()
            results = r.get()  # [sub_obs,sub_obs...]
            obs = pd.concat(results)
            self.adata.obs = obs

            # prefixing
            self._prefixing(col='raw')

        else:
            # compute shaley value
            score_colname = copy.deepcopy(self.total_metrics)
            score_colname.remove('doublet')
            data = np.empty([len(self.query),self.adata.obs.shape[0],len(score_colname)])  # store the metric data for each cell
            '''
            data:
            depth is how many sets of annotations
            height is how many cells
            width is how many score metrics
            '''
            for i,key in enumerate(self.query):
                practical_colname = [name + '@' + key for name in score_colname]
                data[i,:,:] = self.adata.obs[practical_colname].values
            final = []
            intermediate = []

            # computing
            obs = self.adata.obs
            collect = run_shapley(obs,self.query,self.reference,self.size_dict,data,mode,bonus)
            final.extend(collect[0])
            intermediate.extend(collect[1])
            self.adata.obs['final_annotation'] = final
            decisions = list(zip(*intermediate))
            for i,d in enumerate(decisions):
                self.adata.obs['{}_shapley'.format(self.query[i])] = d


            # get raw sctriangulate result
            obs = self.adata.obs
            obs = run_assign(obs)
            self.adata.obs = obs

            # prefixing
            self._prefixing(col='raw')




    def _prefixing(self,col):
        col1 = self.adata.obs[col]
        col2 = self.adata.obs[self.reference]
        col = []
        for i in range(len(col1)):
            concat = self.reference + '@' + col2[i] + '|' + col1[i]
            col.append(concat)
        self.adata.obs['prefixed'] = col


    def pruning(self,method='reassign',discard=None,scale_sccaf=True,layer=None,abs_thresh=10,remove1=True,reference=None,parallel=True,assess_raw=False):
        '''
        Main function. After running ``compute_shapley``, we get **raw** cluster results. Althought the raw cluster is informative,
        there maybe some weired clusters that accidentally win out which doesn't attribute to its biological stability. For example,
        a cluster that only has 3 cells, or very unstable cluster. To ensure the best results, we apply a post-hoc assessment onto the
        raw cluster result, by applying the same set of metrics function to assess the robustness/stability of the raw clusters itself.
        And we will based on that to perform some pruning to get rid of unstable clusters. Finally, the cells within these clusters will be 
        reassigned to thier nearest neighbors.

        :param method: string, valid value: 'reassign', 'rank'. ``rank`` will compute the metrics on all the raw clusters, together with the 
                       ``discard`` parameter which automatically discard clusters ranked at the bottom to remove unstable clusters. ``reassign``
                       will just remove clusters that either has less than ``abs_thresh`` cells or are in the self.invalid attribute list.
        :param discard: int. Least {discard} stable clusters to remove. Default: None, means just rank without removing.
        :param scale_sccaf: boolean. whether to scale the expression data. See ``compute_metrics`` for full explanation. Default: True
        :param abs_thresh: int. clusters have less than {abs_thresh} cells will be discarded in ``reassign`` mode.
        :param remove1: boolean. When reassign the cells in the dicarded clutsers, whether to also reassign the cells who are the only one in each 
                        ``reference`` cluster. Default: True
        :param reference: string. which annotation will serve as the reference.
        :param parallel: boolean, whether to perform this step in parallel. (scatter and gather).Default is True
        :param assess_raw: boolean, whether to run the same set of cluster stability metrics on the raw cluster. Default is False

        Examples::

            sctri.pruning(method='pruning',discard=None)  # just assess and rank the raw clusters
            sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='annotation1')  # remove invalid clusters and reassign the cells within

        '''
        if parallel:
            if method == 'reference':
                obs = reference_pruning(self.adata.obs,self.reference,self.size_dict)
                self.adata.obs = obs

            elif method == 'reassign':
                obs, invalid = reassign_pruning(self,abs_thresh=abs_thresh,remove1=remove1,reference=reference)
                self.adata.obs = obs
                self.invalid = invalid

            elif method == 'rank':
                obs, df = rank_pruning(self,discard=discard,scale_sccaf=scale_sccaf,layer=layer,assess_raw=assess_raw)
                self.adata.obs = obs
                self.uns['raw_cluster_goodness'] = df
                self.adata.obs['confidence'] = self.adata.obs['pruned'].map(df['win_fraction'].to_dict())

        self._prefixing(col='pruned')

        # # finally, generate a celltype sheet
        # obs = self.adata.obs
        # with open(os.path.join(self.dir,'celltype.txt'),'w') as f:
        #     f.write('reference\tcell_cluster\tchoice\n')
        #     for ref,grouped_df in obs.groupby(by=self.reference):
        #         unique = grouped_df['pruned'].unique()
        #         for reassign in unique:
        #             f.write('{}\t{}\n'.format(self.reference + '@' + ref,reassign))



    def get_cluster(self):
        sheet = pd.read_csv(os.path.join(self.dir,'celltype.txt'),sep='\t')
        mapping = {}
        for ref,sub_df in sheet.groupby(by='reference'):
            for cho,subsub_df in sub_df.groupby(by='choice'):
                tmp_list = subsub_df['cell_cluster'].tolist()
                composite_name = ref + '|' + '+'.join(tmp_list)
                for item in tmp_list:
                    original_name = ref + '|' + item
                    mapping[original_name] = composite_name
        self.adata.obs['user_choice'] = self.adata.obs['prefixed'].map(mapping).values


    def elo_rating_like(self):
        '''
        Computing an overall quality score for each annotation, like the idea of ``Elo Rating`` in chess in which it can 
        reflect the probability that one player will win in a chess match. Here, we argue the overall quality score should be defined by 
        the average of all cells' shapley value in all clusters, normalized by the number of players (annotation), and further normalized by
        number of clusters, as our computation of shapley is an additive process such that more player will result in higher shapley value. 
        This step should be run after shapley value were evaluated.

        :return result_dic: a dictionary keyed by annotation, value is overall quality score

        Examples::

            result_dic = sctri.elo_rating_like()

            # {'sctri_rna_leiden_1': 1.5053613872472829, 'sctri_rna_leiden_2': 1.0973714905049967, 'sctri_rna_leiden_3': 1.1032324231884296}
        '''
        obs = self.adata.obs.copy()
        n_p = len(self.query)
        result_dic = {}
        for q in self.query:
            col_cluster_name = q
            col_shapley_name = '_'.join([q,'shapley'])
            elo_rating = 0
            n_c = 0
            for c,sub_df in obs.groupby(by=q):
                elo_rating_c = sub_df[col_shapley_name].mean() / n_p
                elo_rating += elo_rating_c
                n_c += 1
            elo_rating = elo_rating / n_c
            result_dic[q] = elo_rating
        return result_dic
            


        

    def plot_umap(self,col,kind='category',save=True,format='pdf',umap_dot_size=None,umap_cmap='YlOrRd',frameon=False):
        '''
        plotting the umap with either category cluster label or continous metrics are important. Different from the scanpy 
        vanilla plot function, this function automatically generate two umap, one with legend on side and another with legend on data, which
        usually will be very helpful imagine you have > 40 clusters. Secondly, we automatically make all the background dot as light grey, instead of
        dark color.

        :param col: string, which column in self.adata.obs that we want to plot umap from.
        :param kind: string, either 'category' or 'continuous'
        :param save: boolean, whether to save it to disk or not. Default: True
        :param format: string. Which format to save. Default: '.pdf'
        :param umap_dot_size: int/float. the size of dot in scatter plot, if None, using scanpy formula, 120000/n_cells
        :param umap_cmap: string, the matplotlib colormap to use. Default: 'YlOrRd'
        :param frameon: boolean, whether to have the frame on the umap. Default: False

        Examples::

            sctri.plot_umap(col='pruned',kind='category')
            sctri.plot_umap(col='confidence',kind='continous')

        .. image:: ./_static/umap_category.png
            :height: 700px
            :width: 400px
            :align: center
            :target: target      

        .. image:: ./_static/umap_continuous.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target     
        '''
        # col means which column in obs to draw umap on
        if umap_dot_size is None:
            dot_size = 120000/self.adata.obs.shape[0]
        else:
            dot_size = umap_dot_size
        if kind == 'category':
            fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
            sc.pl.umap(self.adata,color=col,frameon=frameon,ax=ax[0],size=dot_size)
            sc.pl.umap(self.adata,color=col,frameon=frameon,legend_loc='on data',legend_fontsize=5,ax=ax[1],size=dot_size)
            if save:
                plt.savefig(os.path.join(self.dir,'umap_sctriangulate_{}.{}'.format(col,format)),bbox_inches='tight')
                plt.close()
        elif kind == 'continuous':
            sc.pl.umap(self.adata,color=col,frameon=frameon,cmap=bg_greyed_cmap(umap_cmap),vmin=1e-5,size=dot_size)
            if save:
                plt.savefig(os.path.join(self.dir,'umap_sctriangulate_{}.{}'.format(col,format)),bbox_inches='tight')
                plt.close()


    def plot_concordance(self,key1,key2,style='3dbar',save=True,format='pdf',cmap=retrieve_pretty_cmap('scphere'),**kwargs):
        '''
        given two annotation, we want to know how the cluster labels from one correponds to the other.

        :param key1: string, first annotation key
        :param key2: string, second annotation key 
        :param style: string, which style of plot, either 'heatmap' or '3dbar'
        :param save: boolean, save the figure or not
        :param format: string, format to save
        :param cmap: string or cmap object, the cmap to use for heatmap

        :return: dataframe, the confusion matrix

        Examples::

            sctri.plot_concordance(key1='azimuth',key2='pruned',style='3dbar')
        
        .. image:: ./_static/3dbar.png
            :height: 350px
            :width: 400px
            :align: center
            :target: target          

        '''
        # construct key1 and key2 map
        key1_map = self.adata.obs[key1].to_frame().groupby(by=key1).apply(lambda x:x.index.to_list()).to_dict()
        key2_map = self.adata.obs[key2].to_frame().groupby(by=key2).apply(lambda x:x.index.to_list()).to_dict()
        # build confusion_df
        confusion_mat = np.empty(shape=(len(key1_map),len(key2_map)),dtype=np.int64)
        for i,(k1,v1) in enumerate(key1_map.items()):
            for j,(k2,v2) in enumerate(key2_map.items()):
                overlap = set(key1_map[k1]).intersection(set(key2_map[k2]))
                n_overlap = len(overlap)
                confusion_mat[i,j] = n_overlap
        confusion_df = pd.DataFrame(data=confusion_mat,index=list(key1_map.keys()),columns=list(key2_map.keys()))
        # plot heatmap
        if style == 'heatmap':
            sns.heatmap(confusion_df,cmap=cmap,**kwargs)  
            if save:
                plt.savefig(os.path.join(self.dir,'concordance_heatmap_{}_{}.{}'.format(key1,key2,format)),bbox_inches='tight')
                plt.close()
        # plot 3D barplot
        elif style == '3dbar':
            fig = plt.figure()
            ax1 = fig.add_subplot(111, projection='3d')
            _x = np.arange(confusion_df.shape[0])
            _y = np.arange(confusion_df.shape[1])
            _xx,_yy = np.meshgrid(_x,_y)
            x,y = _xx.flatten(), _yy.flatten()
            dz = confusion_mat.T.flatten().astype(dtype=np.float64)
            z = np.zeros(len(x))
            dx = np.full(len(x),fill_value=0.6)
            dy = np.full(len(x),fill_value=0.6)
            from scipy.interpolate import interp1d
            from matplotlib import cm
            m = interp1d((dz.min(),dz.max()),(0,255))
            c = [cm.jet(round(i)) for i in m(dz)]
            ax1.bar3d(x, y, z,dx,dy,dz,color=c)
            ax1.set_xlabel('{} cluster labels'.format(key1))
            ax1.set_ylabel('{} cluster labels'.format(key2))
            ax1.set_zlabel('# cells')
            ax1.set_xticks(_x)
            ax1.set_xticklabels(confusion_df.index,fontsize=2)
            ax1.set_yticks(_y+1)
            ax1.set_yticklabels(confusion_df.columns,fontsize=2)
            # # make the panes transparent
            # ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # make the grid lines transparent
            # ax1.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            # ax1.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            # ax1.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
            if save:
                plt.savefig(os.path.join(self.dir,'concordance_3dbarplot_{}_{}.{}'.format(key1,key2,format)),bbox_inches='tight')
                plt.close()

        return confusion_df

    def plot_confusion(self,name,key,save=True,format='pdf',cmap=retrieve_pretty_cmap('scphere'),labelsize=None,**kwargs):
        '''
        plot the confusion as a heatmap.

        :param name: string, either 'confusion_reassign' or 'confusion_sccaf'.
        :param key: string, a annotation name which we want to assess the confusion matrix of the clusters.
        :param save: boolean, whether to save the figure. Default: True.
        :param format: boolean, file format to save. Default: '.pdf'.
        :param cmap: colormap object, Default: scphere_cmap, which defined in colors module.
        :param labelsize: float, this can adjust the label size on yaxis and xaxis for the resultant heatmap
        :param kwargs: additional keyword arguments to sns.heatmap().

        Examples::

            sctri.plot_confusion(name='confusion_reassign',key='sctri_rna_leiden_1')

        .. image:: ./_static/pc.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target          
        '''
        df = self.uns[name][key]
        df = df.apply(func=lambda x:x/x.sum(),axis=1)
        fig,ax = plt.subplots()
        sns.heatmap(df,ax=ax,cmap=cmap,**kwargs) 
        if labelsize is not None: 
            ax.tick_params(labelsize=1)
        if save:
            plt.savefig(os.path.join(self.dir,'confusion_{}_{}.{}'.format(name,key,format)),bbox_inches='tight')
            plt.close()
    
    def plot_cluster_feature(self,key,cluster,feature,enrichment_type='enrichr',save=True,format='pdf'):
        '''
        plot the feature of each single clusters, including:
        
        1. enrichment of artifact genes
        2. marker genes umap
        3. exclusive genes umap
        4. location of clutser umap

        :param key: string. Name of the annation.
        :param cluster: string. Name of the cluster in the annotation.
        :param feature: string, valid value: 'enrichment','marker_genes', 'exclusive_genes', 'location'
        :param enrichmen_type: string, either 'enrichr' or 'gsea'.
        :param save: boolean, whether to save the figure.
        :param format: string, which format for the saved figure.

        Example::

            sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='enrichment')

        .. image:: ./_static/enrichment.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target  

        Example::

            sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='marker_genes')

        .. image:: ./_static/marker_genes.png
            :height: 250px
            :width: 800px
            :align: center
            :target: target  

        Example::

            sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='location')

        .. image:: ./_static/location.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target                     

     '''
        if feature == 'enrichment':
            fig,ax = plt.subplots()
            a = self.uns['marker_genes'][key].loc[cluster,:][enrichment_type]
            ax.barh(y=np.arange(len(a)),width=[item for item in a.values()],color='#FF9A91')
            ax.set_yticks(np.arange(len(a)))
            ax.set_yticklabels([item for item in a.keys()])
            ax.set_title('Marker gene enrichment')
            ax.set_xlabel('-Log10(adjusted_pval)')
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_enrichment.{2}'.format(key,cluster,format)),bbox_inches='tight')
                plt.close()
        elif feature == 'marker_genes':
            a = self.uns['marker_genes'][key].loc[cluster,:]['purify']
            top = a[:10]
            # change cmap a bit
            sc.pl.umap(self.adata,color=top,ncols=5,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_marker_umap.{2}'.format(key,cluster,format)),bbox_inches='tight')
                plt.close()
        elif feature == 'exclusive_genes':
            a = self.uns['exclusive_genes'][key][cluster]  # self.uns['exclusive_genes'][key] is a pd.Series
            a = list(a.keys())
            top = a[:10]
            sc.pl.umap(self.adata,color=top,ncols=5,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_exclusive_umap.{2}'.format(key,cluster,format)),bbox_inches='tight')
                plt.close()
        elif feature == 'location':
            col = [1 if item == str(cluster) else 0 for item in self.adata.obs[key]]
            self.adata.obs['tmp_plot'] = col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap('YlOrRd'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_location_umap.{2}'.format(key,cluster,format)),bbox_inches='tight')
                plt.close()

    def plot_heterogeneity(self,key,cluster,style,col='pruned',save=True,format='pdf',genes=None,umap_zoom_out=True,umap_dot_size=None,
                           subset=None,marker_gene_dict=None,jitter=True,rotation=60,single_gene=None,dual_gene=None,multi_gene=None,merge=None,
                           to_sinto=False,to_samtools=False,cmap='YlOrRd',heatmap_cmap='viridis',heatmap_scale=None,heatmap_regex=None,heatmap_direction='include',
                           heatmap_n_genes=None,heatmap_cbar_scale=None,gene1=None,gene2=None,kind=None,hist2d_bins=50,hist2d_cmap=bg_greyed_cmap('viridis'),
                           hist2d_vmin=1e-5,hist2d_vmax=None,scatter_dot_color='blue',contour_cmap='viridis',contour_levels=None,contour_scatter=True,
                           contour_scatter_dot_size=5,contour_train_kde='valid',surface3d_cmap='coolwarm',**kwarg): 
        '''
        Core plotting function in scTriangulate.

        :param key: string, the name of the annotation.
        :param cluster: string, the name of the cluster.
        :param stype: string, valid values are as below:

            * **umap**: plot the umap of this cluster (including its location and its suggestive heterogeneity)
            * **heatmap**: plot the heatmap of the differentially expressed features across all sub-populations within this cluster.
            * **build**: plot both the umap and heatmap, benefit is the column and raw colorbar of the heatmap is consistent with the umap color
            * **heatmap_custom_gene**, plot the heatmap, but with user-defined gene dictionary.
            * **heatmap+umap**, it is the umap + heatmap_custom_gene, and the colorbars are matching
            * **violin**: plot the violin plot of the specified genes across sub populations.
            * **single_gene**: plot the gradient of single genes across the cluster.
            * **dual_gene**: plot the dual-gene plot of two genes across the cluster, usually these two genes should correspond to the marker genes in two of the sub-populations.
            * **multi_gene**: plot the multi-gene plot of multiple genes across the cluster.
            * **cellxgene**: output the h5ad object which are readily transferrable to cellxgene. It also support atac pseudobuld analysis with ``to_sinto`` or ``to_samtools`` arguments.
            * **sankey**: plot the sankey plot showing fraction/percentage of cells that flow into each sub population, requiring plotly if you only need html, and kaleido if you need static plot, otherwise, a less pretty matplotlib sankey will be plotted
            * **coexpression**: visualize the coexpression pattern of two features, using contour plot or hist2d

        :param col: string, either 'raw' or 'pruned'.
        :param save: boolean, whether to save or not.
        :param foramt: string, which format to save.
        :param genes: list, for violin plot.
        :param umap_zoom_out: boolean, for the umap, whether to zoom out meaning the scale is the same of the whole umap. Zoom in means an amplified version of this cluster.
        :param umap_dot_size: int/float, for the umap.
        :param subset: list, the sub populations we want to keep for plotting.
        :param marker_gene_dict: dict. The custom genes we want the heatmap to display.
        :param jitter: float, for the violin plot.
        :param rotation: int/float, for the violin plot. rotation of the text. Default: 60
        :param single_gene: string, the gene name for single gene plot
        :param dual_gene: list, the dual genes for dual gene plot.
        :param multi_genes: list, the multiple genes for multi gene plot.
        :param merge: nested list, the sub-populations that we want to merge. [('sub-c1','sub-c2'),('sub-c3','sub-c4')]
        :param to_sinto: boolean, for cellxgene mode, output the txt files for running sinto to generate pseudobulk bam files.
        :param to_samtools: boolean,for cellxgene mode, output the txt files for running samtools to generate the pseudobulk bam files.
        :param cmap: a valid string for matplotlib cmap or scTriangulate color module retrieve_pretty_cmap function return object, default is 'YlOrRd', will be used for umap

        The following will be used for heatmap only:

        :param heatmap_cmap: a valid string for matplotlib cmap or scTriangulate color module retrieve_pretty_cmap function return object, default is 'viridis'.
        :param heatmap_scale: None, minmax, median, mean, z_score, default is None, useful when very large or small values exist in the adata.X, scaling can yield better visual effects

            * ``None`` means no scale will be performed, the raw valus shown in adata.X will be plotted in the heatmap
            * ``minmax`` means the raw values will be row-scaled to [0,1] using a MinMaxScaler
            * ``median`` means the raw values will be row-scaled via substracting by the median per row
            * ``mean`` means the raw values will be row-scaled via substracting by the mean per row
            * ``z_score`` means the raw values will be row-scaled via Scaling (mean-centered and variance normalized)

        :param heatmap_regex: None or a raw string for example r'^AB_' (meaning selecing all ADT features as scTriangulate by default prefix ADT features will AB_), the usage of that is to only display certain features from certain modlaities. The underlying implementation is just a regular expression selection. 
        :param heatmap_direction: string, 'include' or 'exclude', it is along with the heatmap_regex parameter, include means doing positive selection, exclude means to exclude the features that match with the heatmap_regex
        :param heatmap_n_genes: an integer, by default, program display 50//n_cluster genes for each cluster, this will overwrite the default.
        :param heatmap_cbar_scale: None or a tuple or a fraction. A tuple for example (-0.5,0.5) will clip the colorbar within -0.5 to 0.5, a fraction number for instance 0.25, will shrink the default colorbar range say -1 to 1 to -0.25 to 0.25

        The following will be used for coexpression plot only:

        :param gene1: the first gene/features to inspect, the gene name, a string.
        :param gene2: the second gene/features to inspect, the gene name, a string.
        :param kind: a string, 'scatter' or 'hist2d' or 'contour' or 'contourf' or 'surface3d', those are all supported figure types to represent the coexpression pattern of two features.
        :param hist2d_bins: integer, default is 50, only used is the kind is hist2d, it will determine the number of the bins
        :param hist2d_cmap: a valid matplotlib cmap string, default is bg_greyed_cmap('viridis'), only used for hist2d
        :param hist2d_vmin: the min value for hist2d graph, default is 1e-5, useful if you want to make the low expressin region lightgrey.
        :param hist2d_vmax: the max value for hist2d graph, default is None
        :param scatter_dot_color: the color of the scatter plot dot, default is 'blue'
        :param contour_cmap: the valid matplotlib camp string, default is 'viridis'
        :param contour_levels: an integer, the levels of contours to show, default is None
        :param contour_scatter: boolean and default is True, whether or not to show the scatter plot on top of the contour plot
        :param contour_scatter_dot_size: float or integer, the dot size of the scatter plot on top of contour plot, default is 5.
        :param contour_train_kde: a string, either 'valid', 'semi-vaid' or 'full', it determines what subset of dots will be used for inferring the kernel

            * ``valid``: only data points that are non-zero for both gene1 and gene2
            * ``semi_valid``: only data points that are non-zero for at least one of the gene
            * ``full``: all data points will be used for kde estimation

        :param surface_3d_cmap: a valid matplotlib cmap string, for surface 3d plot, the default would be 'coolwarm'


        Example::
        
            sctri.plot_heterogeneity('leiden1','0','umap',subset=['leiden1@0','leiden3@10'])
            sctri.plot_heterogeneity('leiden1','0','heatmap',subset=['leiden1@0','leiden3@10'])
            sctri.plot_heterogeneity('leiden1','0','violin',subset=['leiden1@0','leiden3@10'],genes=['MAPK14','ANXA1'])
            sctri.plot_heterogeneity('leiden1','0','sankey')
            sctri.plot_heterogeneity('leiden1','0','cellxgene')
            sctri.plot_heterogeneity('leiden1','0','heatmap+umap',subset=['leiden1@0','leiden3@10'],marker_gene_dict=marker_gene_dict)
            sctri.plot_heterogeneity('leiden1','0','dual_gene',dual_gene=['MAPK14','CD52'])
            sctri.plot_heterogeneity('leiden1','0','coexpression',gene1='MAPK14',genes='CD52',kind='contour')
   
        '''
        adata_s = self.adata[self.adata.obs[key]==cluster,:].copy()
        # remove prior color stamps
        tmp = adata_s.uns
        tmp.pop('{}_colors'.format(col),None)
        adata_s.uns = tmp

        # only consider the sub-populations in subset list
        if subset is not None:
            adata_s = adata_s[adata_s.obs[col].isin(subset),:].copy()
        if merge is not None:
            # if merge is not None, merge the sub-populations that are in each list
            # and make sure it execucate after subetting, so don't contain sub-populations that not in subset.
            # merge argument should be a nested list [('leiden1@3','leiden2@3'),('leiden3@4','leiden4@5')]
            the_map = {}
            # first put all sub_pop that needs to be concated in the map
            for need_merge in merge:
                new_concat_name = '+'.join(need_merge)
                for sub_pop in need_merge:
                    the_map[sub_pop] = new_concat_name
            # then check the remaining pop that doesn't neee to be concated, put into the_map
            all_pop = adata_s.obs[col].unique()
            remain_pop = [item for item in all_pop if item not in the_map.keys()]
            for item in remain_pop:
                the_map[item] = item
            # now map and get new column, and modifiy it back to "col"
            tmp_new_col = adata_s.obs[col].map(the_map).values
            adata_s.obs[col] = tmp_new_col

        if style == 'build':  # draw umap and heatmap
            

            # umap
            fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0.5},figsize=(5,10))
            # ax1
            sc.pl.umap(adata_s,color=[col],ax=axes[0])
            # ax2
            tmp_col = [1 if item == str(cluster) else 0 for item in self.adata.obs[key]]
            self.adata.obs['tmp_plot'] = tmp_col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap(cmap),vmin=1e-5,ax=axes[1])
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,'umap',format)),bbox_inches='tight')
                plt.close()
            self.adata.obs.drop(columns=['tmp_plot'])

            # heatmap
            tmp = adata_s.uns
            tmp.pop('rank_genes_groups',None)
            adata_s.uns = tmp

            if heatmap_scale is not None :   # rowwise scaling in case features from multiple modalities are in differen scale
                if heatmap_scale == 'minmax':
                    from sklearn.preprocessing import MinMaxScaler
                    scaled_X = MinMaxScaler().fit_transform(make_sure_mat_dense(adata_s.X))
                    adata_s.X = scaled_X
                elif heatmap_scale == 'median':
                    scaled_X = make_sure_mat_dense(adata_s.X) - np.median(make_sure_mat_dense(adata_s.X),axis=0)[np.newaxis,:]
                    adata_s.X = scaled_X
                elif heatmap_scale == 'mean':
                    scaled_X = make_sure_mat_dense(adata_s.X) - np.mean(make_sure_mat_dense(adata_s.X),axis=0)[np.newaxis,:]
                    adata_s.X = scaled_X     
                elif heatmap_scale == 'z_score':
                    from sklearn.preprocessing import scale
                    scaled_X = scale(X=make_sure_mat_dense(adata_s.X),axis=0)
                    adata_s.X = scaled_X   

            if len(adata_s.obs[col].unique()) == 1: # it is already unique
                logger_sctriangulate.info('{0} entirely being assigned to one type, no need to do DE'.format(cluster))
                return None
            else:
                sc.tl.rank_genes_groups(adata_s,groupby=col)
                adata_s = filter_DE_genes(adata_s,self.species,self.criterion,heatmap_regex,heatmap_direction)
                number_of_groups = len(adata_s.obs[col].unique())
                if heatmap_n_genes is None:
                    genes_to_pick = 50 // number_of_groups
                else:
                    genes_to_pick = heatmap_n_genes
                if heatmap_cbar_scale is None:   # let scanpy default norm figure that out for you, seems the max and min are not the same as the max/min from the data
                    sc.pl.rank_genes_groups_heatmap(adata_s,n_genes=genes_to_pick,swap_axes=True,key='rank_genes_groups_filtered',cmap=heatmap_cmap)
                else:
                    if isinstance(heatmap_cbar_scale,tuple):
                        v = make_sure_mat_dense(adata_s.X)
                        min_now = heatmap_cbar_scale[0]
                        max_now = heatmap_cbar_scale[1]
                    else:
                        v = make_sure_mat_dense(adata_s.X)
                        max_v = v.max()
                        min_v = v.min()
                        max_v = max([max_v,abs(min_v)])     # make them symmetrical 
                        min_v = max_v * (-1)   
                        max_now = max_v * heatmap_cbar_scale
                        min_now = min_v * heatmap_cbar_scale
                    adata_s.layers['to_plot'] = v                   # very weired fix, have to set a new layer....    
                    sc.pl.rank_genes_groups_heatmap(adata_s,layer='to_plot',n_genes=genes_to_pick,swap_axes=True,key='rank_genes_groups_filtered',cmap=heatmap_cmap,
                                                    vmin=min_now,vmax=max_now)
                    
                if save:
                    plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,'heatmap',format)),bbox_inches='tight')
                    plt.close()

        elif style == 'single_gene':
            fig,ax = plt.subplots()
            if umap_dot_size is None:
                s = 120000/self.adata.obs.shape[0]
            else:
                s = umap_dot_size
            if umap_zoom_out:
                umap_whole = self.adata.obsm['X_umap']
                umap_x_lim = (umap_whole[:,0].min(),umap_whole[:,0].max())
                umap_y_lim = (umap_whole[:,1].min(),umap_whole[:,1].max())
                ax.set_xlim(umap_x_lim)
                ax.set_ylim(umap_y_lim)
            sc.pl.umap(adata_s,color=[single_gene],size=s,ax=ax,cmap=bg_greyed_cmap(cmap),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}_{}.{}'.format(key,cluster,col,style,single_gene,format)),bbox_inches='tight')
                plt.close()


        
        elif style == 'dual_gene':
            if umap_dot_size is None:
                s = 120000/self.adata.obs.shape[0]
            else:
                s = umap_dot_size
            umap_whole = self.adata.obsm['X_umap']
            umap_x_lim = (umap_whole[:,0].min(),umap_whole[:,0].max())
            umap_y_lim = (umap_whole[:,1].min(),umap_whole[:,1].max())
            dual_gene_plot(adata_s,dual_gene[0],dual_gene[1],s=s,save=save,format=format,dir=self.dir,umap_lim=[umap_x_lim,umap_y_lim])

        elif style == 'multi_gene':
            if umap_dot_size is None:
                s = 120000/self.adata.obs.shape[0]
            else:
                s = umap_dot_size            
            umap_whole = self.adata.obsm['X_umap']
            umap_x_lim = (umap_whole[:,0].min(),umap_whole[:,0].max())
            umap_y_lim = (umap_whole[:,1].min(),umap_whole[:,1].max())
            multi_gene_plot(adata_s,multi_gene,s=s,save=save,format=format,dir=self.dir,umap_lim=[umap_x_lim,umap_y_lim])


        elif style == 'heatmap+umap':
            '''first draw umap'''
            fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0.5},figsize=(5,10))
            # ax1
            if umap_zoom_out:
                umap_whole = self.adata.obsm['X_umap']
                umap_x_lim = (umap_whole[:,0].min(),umap_whole[:,0].max())
                umap_y_lim = (umap_whole[:,1].min(),umap_whole[:,1].max())
                axes[0].set_xlim(umap_x_lim)
                axes[0].set_ylim(umap_y_lim)
            if umap_dot_size is None:
                sc.pl.umap(adata_s,color=[col],ax=axes[0],size=120000/self.adata.obs.shape[0])
            else:
                sc.pl.umap(adata_s,color=[col],ax=axes[0],size=umap_dot_size)
            
            # ax2
            if subset is None:
                tmp_col = [1 if item == str(cluster) else 0 for item in self.adata.obs[key]]
            else:
                tmp_col = []
                for i in range(self.adata.obs.shape[0]):
                    ori_cluster_label = self.adata.obs[key][i]
                    prune_cluster_label = self.adata.obs[col][i]
                    if ori_cluster_label == str(cluster) and prune_cluster_label in subset:
                        tmp_col.append(1)
                    else:
                        tmp_col.append(0)
            self.adata.obs['tmp_plot'] = tmp_col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap(cmap),vmin=1e-5,ax=axes[1])
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,'umap',format)),bbox_inches='tight')
                plt.close()
            self.adata.obs.drop(columns=['tmp_plot']) 

            '''then draw heatmap'''  
            sc.pl.heatmap(adata_s,marker_gene_dict,groupby=col,swap_axes=True,dendrogram=True)
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,'heatmap_custom',format)),bbox_inches='tight')
                plt.close()


        elif style == 'umap':
            fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0.5},figsize=(5,10))
            # ax1
            if umap_zoom_out:
                umap_whole = self.adata.obsm['X_umap']
                umap_x_lim = (umap_whole[:,0].min(),umap_whole[:,0].max())
                umap_y_lim = (umap_whole[:,1].min(),umap_whole[:,1].max())
                axes[0].set_xlim(umap_x_lim)
                axes[0].set_ylim(umap_y_lim)
            if umap_dot_size is None:
                sc.pl.umap(adata_s,color=[col],ax=axes[0],size=120000/self.adata.obs.shape[0])
            else:
                sc.pl.umap(adata_s,color=[col],ax=axes[0],size=umap_dot_size)
            
            # ax2
            if subset is None:
                tmp_col = [1 if item == str(cluster) else 0 for item in self.adata.obs[key]]
            else:
                tmp_col = []
                for i in range(self.adata.obs.shape[0]):
                    ori_cluster_label = self.adata.obs[key][i]
                    prune_cluster_label = self.adata.obs[col][i]
                    if ori_cluster_label == str(cluster) and prune_cluster_label in subset:
                        tmp_col.append(1)
                    else:
                        tmp_col.append(0)
            self.adata.obs['tmp_plot'] = tmp_col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap(cmap),vmin=1e-5,ax=axes[1])
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,style,format)),bbox_inches='tight')
                plt.close()
            self.adata.obs.drop(columns=['tmp_plot'])

        elif style == 'heatmap':
            tmp = adata_s.uns
            tmp.pop('rank_genes_groups',None)
            adata_s.uns = tmp

            if heatmap_scale is not None :   # rowwise scaling in case features from multiple modalities are in differen scale
                if heatmap_scale == 'minmax':
                    from sklearn.preprocessing import MinMaxScaler
                    scaled_X = MinMaxScaler().fit_transform(make_sure_mat_dense(adata_s.X))
                    adata_s.X = scaled_X
                elif heatmap_scale == 'median':
                    scaled_X = make_sure_mat_dense(adata_s.X) - np.median(make_sure_mat_dense(adata_s.X),axis=0)[np.newaxis,:]
                    adata_s.X = scaled_X
                elif heatmap_scale == 'mean':
                    scaled_X = make_sure_mat_dense(adata_s.X) - np.mean(make_sure_mat_dense(adata_s.X),axis=0)[np.newaxis,:]
                    adata_s.X = scaled_X       
                elif heatmap_scale == 'z_score':
                    from sklearn.preprocessing import scale
                    scaled_X = scale(X=make_sure_mat_dense(adata_s.X),axis=0)
                    adata_s.X = scaled_X   

            if len(adata_s.obs[col].unique()) == 1: # it is already unique
                logger_sctriangulate.info('{0} entirely being assigned to one type, no need to do DE'.format(cluster))
                return None
            else:
                sc.tl.rank_genes_groups(adata_s,groupby=col)
                adata_s = filter_DE_genes(adata_s,self.species,self.criterion,heatmap_regex,heatmap_direction)
                number_of_groups = len(adata_s.obs[col].unique())
                if heatmap_n_genes is None:
                    genes_to_pick = 50 // number_of_groups
                else:
                    genes_to_pick = heatmap_n_genes
                if heatmap_cbar_scale is None:   # let scanpy default norm figure that out for you, seems the max and min are not the same as the max/min from the data
                    sc.pl.rank_genes_groups_heatmap(adata_s,n_genes=genes_to_pick,swap_axes=True,key='rank_genes_groups_filtered',cmap=heatmap_cmap)
                else:
                    if isinstance(heatmap_cbar_scale,tuple):
                        v = make_sure_mat_dense(adata_s.X)
                        min_now = heatmap_cbar_scale[0]
                        max_now = heatmap_cbar_scale[1]
                    else:
                        v = make_sure_mat_dense(adata_s.X)
                        max_v = v.max()
                        min_v = v.min()
                        max_v = max([max_v,abs(min_v)])     # make them symmetrical 
                        min_v = max_v * (-1)   
                        max_now = max_v * heatmap_cbar_scale
                        min_now = min_v * heatmap_cbar_scale
                    adata_s.layers['to_plot'] = v                   # very weired fix, have to set a new layer....    
                    sc.pl.rank_genes_groups_heatmap(adata_s,layer='to_plot',n_genes=genes_to_pick,swap_axes=True,key='rank_genes_groups_filtered',cmap=heatmap_cmap,
                                                    vmin=min_now,vmax=max_now)
                if save:
                    plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,style,format)),bbox_inches='tight')
                    plt.close()
                # return scanpy marker genes for each sub-populations
                sc_marker_dict = {}  # key is subgroup, value is a df containing markers
                col_dict = {}   # key is a colname, value is a numpy record array
                colnames = ['names','scores','pvals','pvals_adj','logfoldchanges']
                for item in colnames:
                    col_dict[item] = adata_s.uns['rank_genes_groups_filtered'][item]
                for group in adata_s.obs[col].unique():
                    df = pd.DataFrame()
                    for item in colnames:
                        df[item] = col_dict[item][group]
                    df.dropna(axis=0,how='any',inplace=True)
                    df.set_index(keys='names',inplace=True)
                    sc_marker_dict[group] = df
                return sc_marker_dict

        elif style == 'coexpression':
            plot_coexpression(adata_s,gene1=gene1,gene2=gene2,kind=kind,hist2d_bins=hist2d_bins,hist2d_cmap=hist2d_cmap,
                           hist2d_vmin=hist2d_vmin,hist2d_vmax=hist2d_vmax,scatter_dot_color=scatter_dot_color,contour_cmap=contour_cmap,
                           contour_levels=contour_levels,contour_scatter=contour_scatter,contour_scatter_dot_size=contour_scatter_dot_size,
                           contour_train_kde=contour_train_kde,surface3d_cmap=surface3d_cmap,save=True,outdir=self.dir,
                           name='{}_{}_heterogeneity_{}_{}_{}_{}_{}.{}'.format(key,cluster,col,gene1,gene2,style,kind,format))

        elif style == 'heatmap_custom_gene':
            sc.pl.heatmap(adata_s,marker_gene_dict,groupby=col,swap_axes=True,dendrogram=True)
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,style,format)),bbox_inches='tight')
                plt.close()

        elif style == 'violin':
            sc.pl.violin(adata_s,genes,groupby=col,rotation=rotation,jitter=jitter)
            if save:
                genes = '_'.join(genes).replace('/','_')
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}_{}.{}'.format(key,cluster,col,genes,style,format)),bbox_inches='tight')
                plt.close()
                
        elif style == 'cellxgene':
            if save:
                adata_s.write(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.h5ad'.format(key,cluster,col,style)))
            if to_sinto:
                if not os.path.exists(os.path.join(self.dir,'sinto')):
                    os.mkdir(os.path.join(self.dir,'sinto'))
                adata_s.obs[col].to_csv(os.path.join(self.dir,'sinto','{}_{}_heterogeneity_{}_{}_to_sinto_cells.txt'.format(key,cluster,col,style)),sep='\t',header=None)
            if to_samtools:
                if not os.path.exists(os.path.join(self.dir,'samtools')):
                    os.mkdir(os.path.join(self.dir,'samtools'))
                for key_,sub_df in adata_s.obs[col].to_frame().groupby(by=col):
                    sub_df.to_csv(os.path.join(self.dir,'samtools','{}_{}_heterogeneity_{}_{}_to_samtools_{}.txt'.format(key,cluster,col,style,key_)),sep='\t',header=None,columns=[])            
            
            # how to use to_sinto or to_samtools file for visualization in IGV (take bigwig)?
            # 1. if use to_sinto to build pseudobulk
            # <1> make sure you pip install sinto
            # <2> download whole bam file, assume barcode is in CB tag field
            # <3> run the following command:
            #     sinto filterbarcodes -b /path/to/whole_bam.bam \
            #                         -c /sinto/azimuth_CD8_TCM_heterogeneity_pruned_cellxgene_to_sinto_cells.txt \
            #                         -p 30
            # <4> for each bam file, build bam.bai, then run bamCoverage:
            #     bamCoverage -b $1.bam -o $1.bw --normalizeUsing CPM -p max -bs 1 -of bigwig

            # 2. if use to_samtools to build pseudobulk
            # <1> make sure to load samtools/1.13.0
            # <2> download whole bam file, know where the barcode is stored
            # <3> run the following command:
            #     samtools view -@ 30 -b -o subset.bam -D CB:test.txt pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam
            #     samtools index resultant.bam
            #     bamCoverage -b $1.bam -o $1.bw --normalizeUsing CPM -p max -bs 1 -of bigwig            
            
            return adata_s





        elif style == 'sankey':
            try:
                import plotly.graph_objects as go
            except:
                logger_sctriangulate.warning('no plotly or kaleido library, fall back to matplotlib sankey plot')
                # processing the obs
                df = pd.DataFrame()
                df['ref'] = ['ref'+':'+key+'@'+cluster for _ in range(adata_s.obs.shape[0])]   # ref:gs@ERP4
                df['query'] = [item.split('@')[0] for item in adata_s.obs[col]]  # leiden1
                df['cluster'] = [item for item in adata_s.obs[col]]  # leiden1@5
                from matplotlib.sankey import Sankey
                fig,ax = plt.subplots()
                sankey = Sankey(ax=ax,head_angle=120,shoulder=0)
                # gs to query
                info1 = {target:-sub.shape[0]/df.shape[0] for target,sub in df.groupby(by='query')}
                flows1 = [1]
                flows1.extend(list(info1.values()))
                labels1 = [df['ref'].values[0]]
                labels1.extend(list(info1.keys()))
                orientations1 = [0,0]
                orientations1.extend(np.random.choice([-1,1],size=len(info1)-1).tolist())
                print(info1,flows1,labels1,orientations1)
                sankey.add(flows=flows1,labels=labels1,trunklength=4,orientations=orientations1)
                # each query to cluster
                for target,sub in df.groupby(by='query'):
                    prior_index_connect = labels1.index(target)
                    info2 = {cluster3:-subsub.shape[0]/sub.shape[0] for cluster3,subsub in sub.groupby(by='cluster')}
                    flows2 = [-flows1[prior_index_connect]]
                    flows2.extend(list(info2.values()))
                    labels2 = [target]
                    labels2.extend(list(info2.keys()))
                    orientations2 = [0,0]
                    orientations2.extend(np.random.choice([-1,1],size=len(info2)-1).tolist())
                    print(info2,flows2,labels2,orientations2)
                    sankey.add(flows=flows2,labels=labels2,trunklength=4,orientations=orientations2,prior=0,connect=(prior_index_connect,0))
                diagrams = sankey.finish()
                # adjust the text labels
                all_text = []
                for plot in diagrams:
                    all_text.append(plot.text)
                    all_text.extend(plot.texts)
                [item.set_fontsize(2) for item in all_text]
                # from adjustText import adjust_text
                # adjust_text(all_text,arrowprops=dict(arrowstyle='->',color='orange'))
                if save:
                    plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,style,format)),bbox_inches='tight')
                    plt.close()

            else:
                df = pd.DataFrame()
                df['ref'] = ['ref'+':'+key+'@'+cluster for _ in range(adata_s.obs.shape[0])]   # ref:gs@ERP4
                df['query'] = [item.split('@')[0] for item in adata_s.obs[col]]  # leiden1
                df['cluster'] = [item for item in adata_s.obs[col]]  # leiden1@5

                unique_ref = df['ref'].unique().tolist() # not lexicographically sorted, only one
                unique_query = df['query'].unique().tolist()  # not lexicographically sorted
                unique_cluster = df['cluster'].unique().tolist() # not lexicographically sorted

                # get node label and node color
                node_label = unique_ref + unique_query + unique_cluster
                node_color = pick_n_colors(len(node_label))

                # get link information [(source,target,value),(),()]    
                link = []
                for target, sub in df.groupby(by='query'):
                    link_ref2query = (sub['ref'].values[0],target,sub.shape[0])
                    link.append(link_ref2query)
                    for cluster3, subsub in sub.groupby(by='cluster'):
                        link_query2cluster = (target,cluster3,subsub.shape[0])
                        link.append(link_query2cluster)
                link_info = list(zip(*link))
                link_source = [node_label.index(item) for item in link_info[0]]
                link_target = [node_label.index(item) for item in link_info[1]]
                link_value = link_info[2]
                link_color = ['rgba{}'.format(tuple([infer_to_256(item) for item in to_rgb(node_color[i])] + [0.4])) for i in link_source]

            
                # start to draw using plotly and save using kaleido
                node_plotly = dict(pad = 15, thickness = 15,line = dict(color = "black", width = 0.5),label = node_label,color = node_color)
                link_plotly = dict(source=link_source,target=link_target,value=link_value,color=link_color)
                fig = go.Figure(data=[go.Sankey(node = node_plotly,link = link_plotly)])
                fig.update_layout(title_text='{}_{}_heterogeneity_{}_{}'.format(key,cluster,col,style), font_size=6)
                if save:
                    try:
                        fig.write_image(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,style,format)))
                    except:
                        fig.write_html(os.path.join(self.dir,'{}_{}_heterogeneity_{}_{}.{}'.format(key,cluster,col,style,format)),include_plotlyjs='cdn')


    def plot_two_column_sankey(self,left_annotation,right_annotation,opacity=0.6,pad=3,thickness=10,margin=300,text=True,save=True):
        '''
        sankey plot to show the correpondance between two annotation, for example, annotation1 and annotation2, how many cells from each cluster in
        annotation1 will flow to each cluster in annotation2.

        :param left_annotation: a string, the name of the annotation1
        :param right_annotation: a string, the name of the annotation2
        :param opacity: float number, default is 0.6, the opacity of the sankey strips
        :param pad: float number, default is 3, the gap between blocks vertically
        :param thickness: float number, default is 10, the width of each block
        :param margin: the white margin of the sankey plot, large value means the sankey plot will not consume the whole horizontal space (shrinkaged), default is 300
        :param text: whether to show the text or not, default is True, only set to False if you want to have publication quality static figure, because plotly will add a weired background shady effect on the text, not good for publication, so you can fisrt remove text, then add it back youself manually
        :param save: wheter to save or not, default is True.

        Example::

            sctri.plot_two_column_sankey('leiden1','leiden2',margin=5)

        .. image:: ./_static/two_column_sankey.png
            :height: 300px
            :width: 400px
            :align: center
            :target: target

        '''
        import plotly.graph_objects as go
        import kaleido
        df = self.adata.obs.loc[:,[left_annotation,right_annotation]]
        node_label = df[left_annotation].unique().tolist() + df[right_annotation].unique().tolist()
        node_color = pick_n_colors(len(node_label))
        link = []
        for source,sub in df.groupby(by=left_annotation):
            for target,subsub in sub.groupby(by=right_annotation):
                if subsub.shape[0] > 0:
                    link.append((source,target,subsub.shape[0]))
        link_info = list(zip(*link))
        link_source = [node_label.index(item) for item in link_info[0]]
        link_target = [node_label.index(item) for item in link_info[1]]
        link_value = link_info[2]
        link_color = ['rgba{}'.format(tuple([infer_to_256(item) for item in to_rgb(node_color[i])] + [opacity])) for i in link_source]
        node_plotly = dict(pad = pad, thickness = thickness,line = dict(color = "grey", width = 0.1),label = node_label,color = node_color)
        link_plotly = dict(source=link_source,target=link_target,value=link_value,color=link_color)
        if not text:
            fig = go.Figure(data=[go.Sankey(node = node_plotly,link = link_plotly, textfont=dict(color='rgba(0,0,0,0)',size=1))])
        else:
            fig = go.Figure(data=[go.Sankey(node = node_plotly,link = link_plotly)])
        fig.update_layout(title_text='sankey_{}_{}'.format(left_annotation,right_annotation), font_size=6, margin=dict(l=margin,r=margin))
        if save:
            try:
                fig.write_image(os.path.join(self.dir,'two_column_sankey_{}_{}_text_{}.pdf'.format(left_annotation,right_annotation,text))) 
            except:
                fig.write_html(os.path.join(self.dir,'two_column_sankey_{}_{}_text_{}.pdf'.format(left_annotation,right_annotation,text)),include_plotlyjs='cdn')  




    def plot_circular_barplot(self,key,col,save=True,format='pdf'):
        # col can be 'raw' or 'pruned'
        obs = copy.deepcopy(self.adata.obs)
        reference = key
        obs['value'] = np.full(shape=obs.shape[0], fill_value=1)
        obs = obs.loc[:, [reference, col, 'value']]
        print(obs)
        obs4plot = obs.groupby(by=[reference, col])['value'].sum().reset_index()
        print(obs.groupby(by=[reference, col])['value'])
        print(obs.groupby(by=[reference, col])['value'].sum())    
        print(obs.groupby(by=[reference, col])['value'].sum().reset_index())
        cmap = colors_for_set(obs4plot[reference].unique().tolist())
        obs4plot['color'] = obs4plot[reference].map(cmap).values


        # plot layout
        upper_limit = 100
        lower_limit = 30
        outer_label_padding = 4
        inner_label_padding = 2

        # rescale the heights
        maximum = obs4plot['value'].max()
        minimum = obs4plot['value'].min()
        heights = (upper_limit - lower_limit)/(maximum - minimum)*(obs4plot['value'].values-minimum) + lower_limit
        obs4plot['value'] = heights



        # plotting
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, polar=True)
        ax.axis('off')
        width = 2 * np.pi / obs4plot.shape[0]
        angles = [width * (i + 1) for i in np.arange(obs4plot.shape[0])]
        bars = ax.bar(x=angles, height=obs4plot['value'].values, width=width, bottom=lower_limit, linewidth=2,
                    edgecolor='white', color=obs4plot['color'].values)



        # labels
        ax.text(x=0,y=0,s=reference,ha='center',va='center')
        for angle, height, label, ref in zip(angles, obs4plot['value'], obs4plot[col], obs4plot[reference]):
            rotation = np.rad2deg(angle)
            alignment = ''
            if angle >= np.pi/2 and angle < 3*np.pi/2:
                alignment = 'right'
                rotation = rotation + 180
            else:
                alignment = 'left'
            ax.text(x=angle, y=lower_limit + height + outer_label_padding, s=label,ha=alignment,va='center',
                    rotation=rotation, rotation_mode='anchor')  # outer labels
            #ax.text(x=angle, y=lower_limit - inner_label_padding, s=ref, va='center')  # inner labels

        # legend
        import matplotlib.patches as mpatches
        ax.legend(handles=[mpatches.Patch(color=i) for i in cmap.values()], labels=list(cmap.keys()),
                    loc='upper left', bbox_to_anchor=(0, 0), ncol=4, frameon=False, columnspacing=10,
                    title='Reference:{}'.format(reference),borderaxespad=10)

        

        if save:
            plt.savefig(os.path.join(self.dir,'sctri_circular_barplot_{}.{}'.format(col,format)),bbox_inches='tight')
            plt.close()

    def modality_contributions(self,mode='marker_genes',key='pruned',tops=20,regex_dict={'adt':r'^AB_','atac':r'^chr\d{1,2}'}):
        '''
        calculate teh modality contributions for multi modal analysis, the modality contributions of each modality of each cluster means
        the number of features from this modality that made into the top {tops} feature list. Three columns will be added to obs, they are

            * **adt_contribution**
            * **atac_contribution**
            * **rna_contribution**

        :param mode: string, either 'marker_genes' or 'exclusive_genes'.
        :param key: string, any valid categorical column in self.adata.obs
        :param tops: int, the top n features to consider for each cluster.
        :param regex_adt: raw string, the pattern by which the ADT feature will be defined.
        :param regex_atac: raw string ,the pattern by which the atac feature will be defined.

        Examples::

            sctri.modality_contributions(mode='marker_genes',key='pruned',tops=20)

        '''
        # based on how many features make into top list to measure its contribution
        # need to choose the persepctive, default is pruned column
        # will build a three maps (ADT, ATAC, RNA), each of them {c1:0.3,..} 
        # only within modality comparison makes sense
        
        # step1: build several maps
        meta_map = {}
        meta_map['rna'] = {}
        for k,v in regex_dict.items():
            meta_map[k] = {}
        # step2: get features and importance, and see which category each feature belongs to
        for cluster in self.adata.obs[key].unique():
            if mode == 'marker_genes':
                features = self.uns[mode][key].loc[cluster]['purify']
                tops_features = features[:tops]
                importance = np.arange(start=tops,stop=0,step=-1)
            elif mode == 'exclusive_genes':
                features = self.uns[mode][key].loc[cluster]  # a dict
                tops_features = list(features.keys())[:tops]
                importance = list(features.values())[:tops]
            for f,i in zip(tops_features,importance):
                being_assigned = False
                for k,v in regex_dict.items():
                    if re.search(pattern=v,string=f):
                        being_assigned = True
                        try:
                            meta_map[k][cluster] += i
                        except KeyError:
                            meta_map[k][cluster] = 0
                            meta_map[k][cluster] += i
                if not being_assigned:
                    try:
                        meta_map['rna'][cluster] += i
                    except:
                        meta_map['rna'][cluster] = 0
                        meta_map['rna'][cluster] += i
        # step3: write to the obs for modality contribution        
        self.adata.obs[key] = self.adata.obs[key].astype('O').astype('category')
        for k in meta_map.keys():
            self.adata.obs['{}_contribution'.format(k)] = self.adata.obs[key].map(meta_map[k]).fillna(0).astype('float32').values


    def plot_multi_modal_feature_rank(self,cluster,mode='marker_genes',key='pruned',tops=20,
                                    regex_dict={'adt':r'^AB_','atac':r'^chr\d{1,2}'},save=True,format='.pdf'):

        '''
        plot the top features in each clusters, the features are colored by the modality and ranked by the importance.

        :param cluster: string, the name of the cluster.
        :param mode: string, either 'marker_genes' or 'exclusive_genes'
        :param tops: int, top n features to plot.
        :param regex_adt: raw string, the pattern by which the ADT feature will be defined.
        :param regex_atac: raw string ,the pattern by which the atac feature will be defined.
        :param save: boolean, whether to save the figures.
        :param format: string, the format the figure will be saved.

        Examples::

            sctri.plot_multi_modal_feature_rank(cluster='sctri_rna_leiden_2@10')
       
        .. image:: ./_static/plot_multi_modal_feature_rank.png
            :height: 550px
            :width: 500px
            :align: center
            :target: target 
        '''
        if mode == 'marker_genes':
            features = self.uns[mode][key].loc[cluster]['purify']
            tops_features = features[:tops]
            x = np.arange(tops)
            labels = tops_features
            importance = np.arange(start=tops,stop=0,step=-1)
        elif mode == 'exclusive_genes':
            features = self.uns[mode][key].loc[cluster]  # a dict
            tops_features = list(features.keys())[:tops]
            importance = list(features.values())[:tops]
            x = np.arange(tops)
            labels = tops_features
        colors = []
        candidate_colors = pick_n_colors(len(regex_dict)+1)
        for item in labels:
            being_assigned = False
            for i,(k,v) in enumerate(regex_dict.items()):
                if re.search(pattern=v,string=item):
                    being_assigned = True
                    colors.append(candidate_colors[i])
                    break
            if not being_assigned:
                colors.append(candidate_colors[-1])
        fig,ax = plt.subplots()
        ax.bar(x=x,height=importance,width=0.5,color=colors,edgecolor='k')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.tick_params(axis='x',labelsize=6,labelrotation=90)
        ax.set_xlabel('top features')
        ax.set_ylabel('Rank(importance)')
        ax.set_title('{}_{}_{}_{}_features'.format(mode,key,cluster,tops))
        import matplotlib.patches as mpatches
        ax.legend(handles=[mpatches.Patch(color=i) for i in candidate_colors],labels=list(regex_dict.keys())+['rna'],
                    frameon=False,loc='upper left',bbox_to_anchor=(1,1))
        
        if save:
            plt.savefig(os.path.join(self.dir,'sctri_multi_modal_feature_rank_{}_{}_{}_{}.{}'.format(mode,key,cluster,tops,format)),bbox_inches='tight')
            plt.close()


                              
    def plot_multi_modal_feature_fraction(self,cluster,mode='marker_genes',key='pruned',tops=[10,20,30,50],
                                    regex_adt=r'^AB_',regex_atac=r'^chr\d{1,2}',save=True,format='pdf'):
        if mode == 'marker_genes':
            features = self.uns[mode][key].loc[cluster]['purify']
        elif mode == 'exclusive_genes':
            features = self.uns[mode][key].loc[cluster]
        data = {}
        for top in tops:
            top_rna,top_adt,top_atac = 0,0,0
            top_adt_name = []
            top_features = features[:top]
            for item in top_features:
                if re.search(pattern=regex_adt,string=item):
                    top_adt += 1
                    top_adt_name.append(item)
                elif re.search(pattern=regex_atac,string=item):
                    top_atac += 1
                else:
                    top_rna += 1
            assert top_adt + top_atac + top_rna == top
            data[top] = (top_rna,top_adt,top_atac,top_adt_name)
        # plotting
        frac_rna = []
        frac_atac = []
        adt_names = []
        for k,v in data.items():
            frac_rna.append(v[0]/k)
            frac_atac.append(v[2]/k)
            adt_names.append(v[3])
        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(nrows=2, ncols=len(data), height_ratios=(0.3, 0.7), hspace=0,wspace=0)
        axes1 = [fig.add_subplot(gs[0,i]) for i in range(len(data))]
        ax2 = fig.add_subplot(gs[1, :])
        # ax2 is the stacked barplot
        width = 1/(2*len(data))
        ax2.set_xlim([0,1])
        x_coord = [1/(2*len(data)) * (i*2+1) for i in range(len(data))]
        ax2.bar(x_coord,frac_rna,width=width,align='center',bottom=0,label='RNA feature',color='#D56DF2',edgecolor='k')
        ax2.bar(x_coord,frac_atac,width=width,align='center',bottom=frac_rna,label='ATAC feature',color='#3FBF90',edgecolor='k')
        ax2.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
        text_lower = [(item[0]+item[1])/2 for item in zip(np.full(len(data),0),frac_rna)]
        text_upper = [item[0] + 1/2 * item[1] for item in zip(frac_rna,frac_atac)]
        for i in range(len(x_coord)):
            ax2.text(x_coord[i],text_lower[i],'{:.2f}'.format(frac_rna[i]),va='center',ha='center')
            ax2.text(x_coord[i],text_upper[i],'{:.2f}'.format(frac_atac[i]),va='center',ha='center')  
        ax2.set_xticks(x_coord)
        ax2.set_xticklabels(['top{}'.format(str(i)) for i in tops])  
        ax2.set_ylabel('RNA/ATAC fractions')
        # ax1 is the single pie chart in axes1 list
        for i,lis in enumerate(adt_names):
            n = len(lis)
            if n > 0:
                axes1[i].pie(x=[100/n for i in range(n)],labels=lis,frame=True,labeldistance=None)
                axes1[i].axis('equal')
                axes1[i].tick_params(bottom=False,left=False,labelbottom=False,labelleft=False)
            else:
                axes1[i].tick_params(bottom=False,left=False,labelbottom=False,labelleft=False)
        axes1[0].set_ylabel('ADT features')
        axes1[-1].legend(loc='lower right',bbox_to_anchor=(1,1),ncol=len(data),frameon=False)
        fig.suptitle('{}_frac_{}_{}'.format(mode,key,cluster))
        if save:
            stringy_tops = '_'.join([str(item) for item in tops])
            plt.savefig(os.path.join(self.dir,'sctri_multi_modal_feature_frac_{}_{}_{}_{}.{}'.format(mode,key,cluster,stringy_tops,format)),bbox_inches='tight')
            plt.close()


    def plot_long_heatmap(self,clusters=None,key='pruned',n_features=5,mode='marker_genes',cmap='viridis',save=True,format='pdf',figsize=(6,4.8),
                          feature_fontsize=3,cluster_fontsize=5,heatmap_regex=None,heatmap_direction='include'):

        '''
        the default scanpy heatmap is not able to support the display of arbitrary number of marker genes for each clusters, the max feature is 50.
        this heatmap allows you to specify as many marker genes for each cluster as possible, and the gene name will all the displayed.

        :param clusters: list, what clusters we want to consider under a certain annotation.
        :param key: string, annotation name.
        :param n_features: int, the number of features to display.
        :param mode: string, either 'marker_genes' or 'exclusive_genes'.
        :param cmap: string, matplotlib cmap string.
        :param save: boolean, whether to save or not.
        :param format: string, which format to save.
        :param figsize: tuple, the width and the height of the plot.
        :param feature_fontsize: int/float. the fontsize for the feature.
        :param cluster_fontsize: int/float, the fontsize for the cluster.
        :param heatmap_regex: None or a raw string for example r’^AB_’ (meaning selecing all ADT features as scTriangulate by default prefix ADT features will AB_), the usage of that is to only display certain features from certain modlaities. The underlying implementation is just a regular expression selection.
        :param heatmap_direction: string, ‘include’ or ‘exclude’, it is along with the heatmap_regex parameter, include means doing positive selection, exclude means to exclude the features that match with the heatmap_regex

        Examples::

            sctri.plot_long_umap(n_features=20,figsize=(20,20))

        .. image:: ./_static/long_heatmap.png
            :height: 550px
            :width: 550px
            :align: center
            :target: target         

        '''
        df = self.uns[mode][key]
        # if heatmap_regex and heatmap_direction are present, try to filter the marker genes first
        if heatmap_regex is not None:
            if heatmap_direction == 'include':
                new_purify_col = []
                for lis in df['purify']:
                    new_lis = []
                    for item in lis:
                        pat = re.compile(heatmap_regex)
                        if re.search(pat,item):
                            new_lis.append(item)
                    new_purify_col.append(new_lis)
            elif heatmap_direction == 'exclude':
                new_purify_col = []
                for lis in df['purify']:
                    new_lis = []
                    for item in lis:
                        pat = re.compile(heatmap_regex)
                        if not re.search(pat,item):
                            new_lis.append(item)
                    new_purify_col.append(new_lis) 
            df['purify'] = new_purify_col
        # get feature pool
        ignore_clusters = []
        feature_pool = []
        for i in range(df.shape[0]):
            cluster = df.index[i]
            if len(df.iloc[i]['purify']) == 0:
                ignore_clusters.append(cluster)
                print(color_stdout('{} only has {} markers with the regex specified, this cluster will not be plotted'.format(cluster,len(df.iloc[i]['purify'])),'red'))
                continue
            elif n_features > len(df.iloc[i]['purify']):
                features = df.iloc[i]['purify']
                print('{} only has {} markers with the regex specified, only these markers will be plotted'.format(cluster,len(df.iloc[i]['purify'])))
                continue
            else:
                features = df.iloc[i]['purify'][:n_features]
            feature_pool.extend(features)
        # determine cluster order
        if clusters is None:
            clusters = list(set(df.index).difference(set(ignore_clusters)))
        core_adata = self.adata[self.adata.obs[key].isin(clusters),feature_pool]
        core_df = pd.DataFrame(data=make_sure_mat_dense(core_adata.copy().X),
                               index=core_adata.obs_names,
                               columns=core_adata.var_names)
        core_df['label'] = core_adata.obs[key].values
        centroid_df = core_df.groupby(by='label').apply(lambda x:x.iloc[:,:-1].mean(axis=0))
        dense_distance_mat = pdist(centroid_df.values,'euclidean')
        linkage_mat = linkage(dense_distance_mat,method='ward',metric='enclidean')
        leaf_order = leaves_list(linkage_mat)
        cluster_order = [centroid_df.index[i] for i in leaf_order]
        # relationship feature-cluster and barcode-cluster, and vice-versa
        feature_cluster_df = pd.DataFrame({'feature':[],'cluster':[]})
        for i in range(df.shape[0]):
            cluster = df.index[i]
            if cluster in ignore_clusters:
                continue
            if n_features > len(df.iloc[i]['purify']):
                features = df.iloc[i]['purify']
            else:
                features = df.iloc[i]['purify'][:n_features]   
            chunk = pd.DataFrame({'feature':features,'cluster':np.full(len(features),fill_value=cluster)})
            feature_cluster_df = pd.concat([feature_cluster_df,chunk],axis=0) 
        feature_to_cluster = feature_cluster_df.groupby(by='feature')['cluster'].apply(lambda x:x.values[0]).to_dict()
        cluster_to_feature = feature_cluster_df.groupby(by='cluster')['feature'].apply(lambda x:x.tolist()).to_dict()

        barcode_cluster_df = pd.DataFrame({'barcode':core_adata.obs_names.tolist(),'cluster':core_adata.obs[key]})
        barcode_to_cluster = barcode_cluster_df.groupby(by='barcode')['cluster'].apply(lambda x:x.values[0]).to_dict()
        cluster_to_barcode = barcode_cluster_df.groupby(by='cluster')['barcode'].apply(lambda x:x.tolist()).to_dict()
        # plotting
        fig = plt.figure(figsize=figsize)
        gs = mpl.gridspec.GridSpec(nrows=2,ncols=2,width_ratios=(0.97,0.03),height_ratios=(0.97,0.03),wspace=0.02,hspace=0.02)
        ax1 = fig.add_subplot(gs[0,0])  # heatmap
        ax2 = fig.add_subplot(gs[1,0])  # column cell color bars
        ax3 = fig.add_subplot(gs[0,1])  # row feature color bars
        # ax1, heatmap
        p_feature = []
        for c in cluster_order:
            p_feature.extend(cluster_to_feature[c])
        p_cell = []
        for c in cluster_order:
            p_cell.extend(cluster_to_barcode[c])
        p_adata = self.adata[p_cell,p_feature].copy()
        draw_data = make_sure_mat_dense(p_adata.X).T
        im = ax1.imshow(X=draw_data,cmap=cmap,aspect='auto',interpolation='none')
        ax1.set_xticks([])
        ax1.set_yticks(np.arange(draw_data.shape[0]))
        ax1.set_yticklabels(p_adata.var_names.tolist(),fontsize=feature_fontsize)  
        # ax2, column cell color bars
        p_adata.obs['plot_cluster'] = p_adata.obs_names.map(barcode_to_cluster)
        tmp_frac = [np.count_nonzero(p_adata.obs['plot_cluster'].values==c)/p_adata.obs.shape[0] for c in cluster_order]
        tmp_cum = np.cumsum(tmp_frac)
        x_coords = [(tmp_cum[i] - tmp_frac[i]*1/2) * p_adata.obs.shape[0] for i in range(len(cluster_order))]
        anno_to_color = colors_for_set(np.sort(p_adata.obs['plot_cluster'].unique()))
        cell_column_cbar_mat = p_adata.obs['plot_cluster'].map(anno_to_color).values.reshape(1,-1)
        cell_column_cbar_mat_rgb = hex2_to_rgb3(cell_column_cbar_mat)
        ax2.imshow(X=cell_column_cbar_mat_rgb,aspect='auto',interpolation='none')
        ax2.set_xticks(x_coords)
        ax2.set_xticklabels(cluster_order,rotation=90,fontsize=cluster_fontsize)
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        # ax3, row feature color bars
        p_adata.var['plot_cluster'] = p_adata.var_names.map(feature_to_cluster)
        feature_row_cbar_mat = p_adata.var['plot_cluster'].map(anno_to_color).values.reshape(-1,1)
        feature_row_cbar_mat_rgb = hex2_to_rgb3(feature_row_cbar_mat)
        ax3.imshow(X=feature_row_cbar_mat_rgb,aspect='auto',interpolation='none')
        ax3.tick_params(bottom=False,left=False,labelbottom=False,labelleft=False)
        # add white vline
        s,e = ax1.get_xlim()
        vline_coords = tmp_cum * (e-s) + s
        for x in vline_coords:
            ax1.axvline(x,ymin=0,ymax=1,color='white',linewidth=0.01) 
        # colorbar
        gs.update(right=0.8)
        gs_cbar = mpl.gridspec.GridSpec(nrows=1,ncols=1,left=0.85,top=0.3)
        ax4 = fig.add_subplot(gs_cbar[0,0])
        plt.colorbar(im,cax=ax4)
    
        if save:
            plt.savefig(os.path.join(self.dir,'sctri_long_heatmap_{}.pdf'.format(key)),bbox_inches='tight')
            plt.close()
        # return that can be imported to morpheus
        export = pd.DataFrame(data=draw_data,columns=p_adata.obs_names,index=p_adata.var_names)
        return export


    def _atomic_viewer_figure(self,key):
        for cluster in self.adata.obs[key].unique():
            try:
                self.plot_cluster_feature(key,cluster,'enrichment','enrichr',True,'png')
                self.plot_cluster_feature(key,cluster,'marker_genes','enrichr',True,'png')
                self.plot_cluster_feature(key,cluster,'exclusive_genes','enrichr',True,'png')
                self.plot_cluster_feature(key,cluster,'location','enrichr',True,'png')
            except KeyError:  # the cluster only have one cell, so not in adata_compute when calculating metrics
                continue


    def _atomic_viewer_hetero(self,key,format='png',heatmap_scale=False,heatmap_cmap='viridis',heatmap_regex=None,heatmap_direction='include',
                              heatmap_n_genes=None,heatmap_cbar_scale=None):
        for cluster in self.adata.obs[key].unique():
            self.plot_heterogeneity(key,cluster,'build',format=format,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,
                                    heatmap_direction=heatmap_direction,heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)


    def viewer_cluster_feature_figure(self,parallel=False,select_keys=None,other_umap=None):
        '''
        Generate all the figures for setting up the viewer cluster page.

        :param parallel: boolean, whether to run it in parallel, only work in some linux system, so recommend to not set to True.
        :param select_keys: list, what annotations' cluster we want to inspect.
        :param other_umap: ndarray,replace the umap with another set.

        Examples::

            sctri.viewer_cluster_feature_figure(parallel=False,select_keys=['annotation1','annotation2'],other_umap=None)

        '''
        logger_sctriangulate.info('Building viewer requires generating all the necessary figures, may take several minutes')
        # see if needs to change the umap embedding
        if other_umap is not None:
            ori_umap = self.adata.obsm['X_umap']
            self.adata.obsm['X_umap'] = other_umap
        # create a folder to store all the figures
        if not os.path.exists(os.path.join(self.dir,'figure4viewer')):
            os.mkdir(os.path.join(self.dir,'figure4viewer'))
        ori_dir = self.dir
        new_dir = os.path.join(self.dir,'figure4viewer')
        self.dir = new_dir
        # generate all the figures
        '''doublet plot'''
        self.plot_umap('doublet_scores','continuous',True,'png')
        if platform.system() == 'Linux' and parallel:    # can parallelize
            cores1 = mp.cpu_count()
            cores2 = len(self.cluster)
            cores = min(cores1,cores2)
            pool = mp.Pool(processes=cores)
            logger_sctriangulate.info('spawn {} sub processes for viewer cluster feature figure generation'.format(cores))
            raw_results = [pool.apply_async(func=self._atomic_viewer_figure,args=(key,)) for key in self.cluster.keys()]
            pool.close()
            pool.join()
        else:                               # Windows and Darwin can not parallelize if plotting
            if select_keys is None:
                for key in self.cluster.keys():
                    self._atomic_viewer_figure(key)
            else:
                for key in select_keys:
                    self._atomic_viewer_figure(key)
        # dial back the dir and umap
        self.dir = ori_dir 
        if other_umap is not None:
            self.adata.obsm['X_umap'] = ori_umap

    def viewer_cluster_feature_html(self):
        '''
        Setting up the viewer cluster page.

        Examples::

            sctri.viewer_cluster_feature_html()

        '''
        # create a folder to store all the figures
        if not os.path.exists(os.path.join(self.dir,'figure4viewer')):
            os.mkdir(os.path.join(self.dir,'figure4viewer'))
        # generate html
        with open(os.path.join(self.dir,'figure4viewer','viewer.html'),'w') as f:
            f.write(to_html(self.cluster,self.score,self.total_metrics))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/viewer.js'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/viewer.css'),os.path.join(self.dir,'figure4viewer')))

    def viewer_heterogeneity_figure(self,key,other_umap=None,format='png',heatmap_scale=False,heatmap_cmap='viridis',heatmap_regex=None,heatmap_direction='include',
                                    heatmap_n_genes=None,heatmap_cbar_scale=None):
        '''
        Generating the figures for the viewer heterogeneity page

        :param key: string, which annotation to inspect the heterogeneity.
        :param other_umap: ndarray, replace with other umap embedding.

        Examples::

            sctri.viewer_heterogeneity_figure(key='annotation1',other_umap=None)
        '''
        logger_sctriangulate.info('Building viewer requires generating all the necessary figures, may take several minutes')
        # see if needs to change umap embedding
        if other_umap is not None:
            ori_umap = self.adata.obsm['X_umap']
            self.adata.obsm['X_umap'] = other_umap            
        # create a folder to store all the figures
        if not os.path.exists(os.path.join(self.dir,'figure4viewer')):
            os.mkdir(os.path.join(self.dir,'figure4viewer'))
        else:  # if already exsiting figure4viewer folder, need to clean previous figures for the specific key
            os.system('rm {}'.format(os.path.join(self.dir,'figure4viewer','{}_*_heterogeneity_*'.format(key))))          
        ori_dir = self.dir
        new_dir = os.path.join(self.dir,'figure4viewer')
        self.dir = new_dir
        self._atomic_viewer_hetero(key,format=format,heatmap_scale=heatmap_scale,heatmap_cmap=heatmap_cmap,heatmap_regex=heatmap_regex,heatmap_direction=heatmap_direction,
                                   heatmap_n_genes=heatmap_n_genes,heatmap_cbar_scale=heatmap_cbar_scale)
        # dial back
        self.dir = ori_dir
        if other_umap is not None:
            self.adata.obsm['X_umap'] = ori_umap

    def viewer_heterogeneity_html(self,key):
        '''
        Setting up viewer heterogeneity page

        :param key: string, which annotation to inspect.

        Examples::

            sctri.viewer_heterogeneity_html(key='annotation1')

        '''
        # create a folder to store all the figures
        if not os.path.exists(os.path.join(self.dir,'figure4viewer')):
            os.mkdir(os.path.join(self.dir,'figure4viewer'))
        key_cluster_dict = copy.deepcopy(self.cluster)
        if key not in key_cluster_dict.keys():
            key_cluster_dict[key] = self.adata.obs[key].unique().tolist()
        with open(os.path.join(self.dir,'figure4viewer','inspection_{}.html'.format(key)),'w') as f:
            f.write(inspection_html(key_cluster_dict,key)) 
        # first copy      
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/inspection.js'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/inspection.css'),os.path.join(self.dir,'figure4viewer')))






# ancillary functions for main class
def penalize_artifact_void(obs,query,stamps,metrics):
    '''
    penalize_artifact_void core function
    '''
    for stamp in stamps:
        metrics_cols = obs.loc[:,[item2+'@'+item1 for item1 in query for item2 in metrics]]
        cluster_cols = obs.loc[:,query]
        df = cluster_cols.apply(func=lambda x:pd.Series(data=[x.name+'@'+str(item) for item in x],name=x.name),axis=0)
        df_repeat = pd.DataFrame(np.repeat(df.values,len(metrics),axis=1))
        truth = pd.DataFrame(data=(df_repeat == stamp).values,index=metrics_cols.index,columns=metrics_cols.columns)
        tmp = metrics_cols.mask(truth,0)
        obs.loc[:,[item2+'@'+item1 for item1 in query for item2 in metrics]] = tmp
    return obs



def each_key_run(sctri,key,scale_sccaf,layer,added_metrics_kwargs=None):
    folder = sctri.dir
    adata = sctri.adata   # here modify adata will still affect sctri.adata
    species = sctri.species
    criterion = sctri.criterion
    metrics = sctri.metrics
    add_metrics = sctri.add_metrics
    total_metrics = sctri.total_metrics

    try:
        assert issparse(adata.X) == False
    except AssertionError:
        adata.X = adata.X.toarray()  

    # remove cluster that only have 1 cell, for DE analysis, adata_to_compute is just a view of adata
    adata_to_compute = check_filter_single_cluster(adata,key)    # here adata_to_compute is just a view of original adata, how I guarantee adata won't change?
                                                                 # it is in each stability function, because potential modification of the adata_to_compute, I will
                                                                 # make a copy and delete it once the computation is done (garbage collection), so that the view (adata_to_compute)
                                                                 # won't be changed, and the adata and sctri.adata won't be changed either.

    # a dynamically named dict
    cluster_to_metric = {}
    '''marker gene'''
    marker_genes = marker_gene(adata_to_compute,key,species,criterion,folder)
    logger_sctriangulate.info('Process {}, for {}, finished marker genes finding'.format(os.getpid(),key))
    '''reassign score'''
    cluster_to_metric['cluster_to_reassign'], confusion_reassign = reassign_score(adata_to_compute,key,marker_genes)
    logger_sctriangulate.info('Process {}, for {}, finished reassign score computing'.format(os.getpid(),key))
    '''tfidf10 score'''
    cluster_to_metric['cluster_to_tfidf10'], exclusive_genes = tf_idf10_for_cluster(adata_to_compute,key,species,criterion,layer=layer)
    logger_sctriangulate.info('Process {}, for {}, finished tfidf score computing'.format(os.getpid(),key))
    '''SCCAF score'''
    cluster_to_metric['cluster_to_SCCAF'], confusion_sccaf = SCCAF_score(adata_to_compute,key, species, criterion,scale_sccaf)
    logger_sctriangulate.info('Process {}, for {}, finished SCCAF score computing'.format(os.getpid(),key))
    '''doublet score'''
    cluster_to_metric['cluster_to_doublet'] = doublet_compute(adata_to_compute,key)
    logger_sctriangulate.info('Process {}, for {}, finished doublet score assigning'.format(os.getpid(),key))
    '''added other scores'''
    for (metric,func),single_kwargs in zip(add_metrics.items(),added_metrics_kwargs):
        cluster_to_metric['cluster_to_{}'.format(metric)] = func(adata_to_compute,key,**single_kwargs)
        logger_sctriangulate.info('Process {}, for {}, finished {} score computing'.format(os.getpid(),key,metric))


    collect = {'key':key}  # collect will be retured to main program
    '''collect all default metrics and added metrics'''
    for metric in total_metrics:
        collect['col_{}'.format(metric)] = adata.obs[key].astype('str').map(cluster_to_metric['cluster_to_{}'.format(metric)]).fillna(0).values
    '''collect score info and cluster info'''
    score_info = cluster_to_metric  # {cluster_to_reassign:{cluster1:0.45}}
    cluster_info = list(cluster_to_metric['cluster_to_reassign'].keys())  #[cluster1,cluster2,cluster3]
    collect['score_info'] = score_info
    collect['cluster_info'] = cluster_info
    '''collect uns including genes and confusion matrix'''
    collect['marker_genes'] = marker_genes
    collect['exclusive_genes'] = exclusive_genes
    collect['confusion_reassign'] = confusion_reassign
    collect['confusion_sccaf'] = confusion_sccaf

    del adata

    return collect


def run_shapley(obs,query,reference,size_dict,data,mode,bonus):
    logger_sctriangulate.info('process {} need to process {} cells for shapley computing'.format(os.getpid(),data.shape[1]))
    final = []
    intermediate = []
    for i in range(data.shape[1]):
        layer = data[:,i,:]
        result = []
        for j in range(layer.shape[0]):
            result.append(wrapper_shapley(j,layer,mode,bonus))
        cluster_row = obs.iloc[i].loc[query].values
        to_take = which_to_take(result,query,reference,cluster_row,size_dict)   # which annotation this cell should adopt
        final.append(to_take)    
        intermediate.append(result)
    return final,intermediate


def run_assign(obs):  
    logger_sctriangulate.info('process {} need to process {} cells for raw sctriangulte result'.format(os.getpid(),obs.shape[0]))   
    assign = []
    for i in range(obs.shape[0]):
        name = obs.iloc[i,:].loc['final_annotation']
        cluster = obs.iloc[i,:].loc[name]
        concat = name + '@' + cluster
        assign.append(concat)   
    obs['raw'] = assign
    return obs

def filter_DE_genes(adata,species,criterion,regex=None,direction='include'):
    de_gene = pd.DataFrame.from_records(adata.uns['rank_genes_groups']['names']) #column use field name, index is none by default, so incremental int value
    # first filter out based on the level of artifact genes
    artifact = set(read_artifact_genes(species,criterion).index)
    de_gene.mask(de_gene.isin(artifact),inplace=True)
    # second filter based on regex and direction (include or exclude)
    if regex is not None:
        pat = re.compile(regex)
        if direction == 'include':
            de_gene.where(de_gene.applymap(lambda x: True if pd.isna(x) or re.search(pat,x) else False),inplace=True)
        elif direction == 'exclude':
            de_gene.mask(de_gene.applymap(lambda x: False if pd.isna(x) or not re.search(pat,x) else True),inplace=True)
    adata.uns['rank_genes_groups_filtered'] = adata.uns['rank_genes_groups'].copy()
    adata.uns['rank_genes_groups_filtered']['names'] = de_gene.to_records(index=False)
    return adata












