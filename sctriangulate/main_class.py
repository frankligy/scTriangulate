import sys
import os
import copy
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import seaborn as sns
from anytree import Node, RenderTree
from scipy.sparse import issparse,csr_matrix
import multiprocessing as mp
import platform
import logging
import subprocess

import scanpy as sc
import anndata as ad
import gseapy as gp
import scrublet as scr

from .shapley import *
from .metrics import *
from .viewer import *
from .prune import *
from .colors import *






# define ScTriangulate Object
class ScTriangulate(object):

    def __init__(self,dir,adata,query,reference,species='mice',criterion=2,verbose=1):

        self.verbose = verbose
        self.dir = dir
        self._create_dir_if_not_exist()
        self.adata = adata
        self.query = query
        self.reference = reference
        self.species = species
        self.criterion = criterion
        self.score = {}
        self.cluster = {}
        self.uns = {}
        self.metrics = ['reassign','tfidf10','SCCAF','doublet']   # default metrics
        self.add_metrics = {}                           # user can add their own, key is metric name, value is callable
        self.total_metrics = self.metrics               # all metrics considered

        self._set_logging()          
        self._check_adata()
        self.size_dict, _ = get_size(self.adata.obs,self.query)
        self.invalid = []





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
        invalid_chars = ['/','@','$']
        if self.reference in self.query:
            all_keys = self.query
        else:
            all_keys = copy.deepcopy(self.query)
            all_keys.append(self.reference)
        for key in all_keys:
            for ichar in invalid_chars:
                self.adata.obs[key] = self.adata.obs[key].str.replace(ichar,'_')
        
        ## replace key name
        for key in all_keys:
            for ichar in invalid_chars:
                self.adata.obs.rename(columns={key:key.replace(ichar,'_')},inplace=True)   
        # step3: remove index name for smooth h5ad writing
        self.adata.obs.index.name = None
        self.adata.var.index.name = None

    def _set_logging(self):
        # get all logger
        global logger_sctriangulate
        logger_sctriangulate = logging.getLogger(__name__)
        logger_scanpy = logging.getLogger('scanpy')
        logger_gseapy = logging.getLogger('gseapy')
        logger_scrublet = logging.getLogger('scrublet')

        # make other logger silent
        logger_scanpy.setLevel(logging.ERROR)
        logger_gseapy.setLevel(logging.ERROR)
        logger_scrublet.setLevel(logging.ERROR)

        # configure own logger
        if self.verbose == 1:
            c_handler = logging.StreamHandler()
            c_handler.setLevel(logging.INFO)
            c_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s' )
            c_handler.setFormatter(c_formatter)        
            logger_sctriangulate.addHandler(c_handler)
            logger_sctriangulate.info('choosing console logging')

        elif self.verbose == 2:
            f_handler = logging.FileHandler(os.path.join(self.dir,'scTriangulate.log'))
            f_handler.setLevel(logging.INFO)
            f_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s' )
            f_handler.setFormatter(f_formatter)
            logger_sctriangulate.addHandler(f_handler)
            logger_sctriangulate.info('choosing file logging')


    def _to_dense(self):
        self.adata.X = self.adata.X.toarray() 
        
    def _to_sparse(self):
        self.adata.X = csr_matrix(self.adata.X)

    def obs_to_df(self,name='sctri_inspect_obs.txt'):
        self.adata.obs.to_csv(os.path.join(self.dir,name),sep='\t')

    def var_to_df(self,name='sctri_inspect_var.txt'):
        self.adata.var.to_csv(os.path.join(self.dir,name),sep='\t')

    def gene_to_df(self,mode,key):
        '''mode is marker_genes or exclusive_genes'''
        self.uns['{}'.format(mode)][key].to_csv(os.path.join(self.dir,'sctri_gene_to_df_{}_{}.txt'.format(mode,key)),sep='\t')

    def confusion_to_df(self,mode,key):
        '''mode is confusion_reassign or confusion_sccaf'''
        self.uns['{}'.format(mode)][key].to_csv(os.path.join(self.dir,'sctri_confusion_to_df_{}_{}.txt'.format(mode,key)),sep='\t')

    def get_metrics_and_shapley(self,barcode,save=False):
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

    def serialize(self,name='sctri_pickle.p'):
        with open(os.path.join(self.dir,name),'wb') as f:
            pickle.dump(self,f)

    @staticmethod
    def deserialize(name):
        with open(name,'rb') as f:
            sctri = pickle.load(f)
        sctri._set_logging()
        logger_sctriangulate.info('unpickled {} to memory'.format(name))
        return sctri

    def add_new_metrics(self,add_metrics):
        for metric,func in add_metrics.items():
            self.add_metrics[metric] = func
        self.total_metrics.extend(list(self.add_metrics.keys()))

    def winners_statistics(self,col,plot,save):
        new_size_dict = {}  # {gs@ERP4: 100}
        for key,value in self.size_dict.items():
            for sub_key,sub_value in value.items():
                composite_key = key + '@' + sub_key
                composite_value = sub_value
                new_size_dict[composite_key] = composite_value
        obs = self.adata.obs
        winners = obs[col]
        winners_vc = winners.value_counts()
        winners_size = winners_vc.index.to_series().map(new_size_dict)
        winners_prop = winners_vc / winners_size
        winners_stats = pd.concat([winners_vc,winners_size,winners_prop],axis=1)
        winners_stats.columns = ['counts','size','proportion']
        winners_stats.sort_values(by='proportion',inplace=True)
        self.winners_stats = winners_stats
        if plot:
            a = winners_stats['proportion']
            fig,ax = plt.subplots()
            ax.barh(y=np.arange(len(a)),width=[item for item in a.values],color='#FF9A91')
            ax.set_yticks(np.arange(len(a)))
            ax.set_yticklabels([item for item in a.index],fontsize=3)
            ax.set_title('Winners statistics')
            ax.set_xlabel('proportion of clusters that win')
            if save:
                plt.savefig(os.path.join(self.dir,'winners_statistics.pdf'),bbox_inches='tight')
                plt.close()

    def clusterability(self,col,plot,save):
        bucket = {}   # {ERP4:5}
        obs = self.adata.obs
        for ref,grouped_df in obs.groupby(by=self.reference):
            unique = grouped_df[col].unique()
            bucket[ref] = len(unique)
        bucket = {k: v for k, v in sorted(bucket.items(), key=lambda x: x[1])}
        if plot:
            fig,ax = plt.subplots()
            ax.scatter(x=np.arange(len(bucket)),y=list(bucket.values()),c='k',s=15)
            ax.set_xticks(np.arange(len(bucket)))
            ax.set_xticklabels(list(bucket.keys()),fontsize=3,rotation=90)
            ax.set_title('{} clusterablity'.format(self.reference))
            ax.set_ylabel('clusterability: # sub-clusters')
            if save:
                plt.savefig(os.path.join(self.dir,'{}_clusterability.pdf'.format(self.reference)),bbox_inches='tight')
                plt.close()


    def display_hierarchy(self,col):
        obs = self.adata.obs
        root = Node(self.reference)
        hold_ref_var = {}
        for ref,grouped_df in obs.groupby(by=self.reference):
            hold_ref_var[ref] = Node(ref,parent=root)
            unique = grouped_df[col].unique()
            if len(unique) == 1: # no sub-clusters
                continue
            else:
                hold_cluster_var = {}
                for item in unique:
                    hold_cluster_var[item] = Node(item,parent=hold_ref_var[ref])
        with open(os.path.join(self.dir,'display_hierarchy_{}_{}.txt'.format(self.reference,col)),'a') as f:
            for pre, fill, node in RenderTree(root):
                print("%s%s" % (pre, node.name),file=f)


    def prune_statistics(self,print=False):
        obs = self.adata.obs
        raw = obs['raw']
        pruned = obs['pruned']
        raw_vc = raw.value_counts()
        pruned_vc = pruned.value_counts()
        pruned_vc_dict = pruned_vc.to_dict()
        tmp = raw_vc.index.map(pruned_vc_dict).fillna(value=0)
        stats_df = raw_vc.to_frame()
        stats_df['pruned'] = tmp.values
        stats_df.sort_values(by='pruned',inplace=True,ascending=False)
        self.prune_stats = stats_df
        if print:
            self.prune_stats.to_csv(os.path.join(self.dir,'sctri_prune_statistics.txt'),sep='\t')


    def doublet_predict(self):
        if issparse(self.adata.X):
            self._to_dense()
        counts_matrix = self.adata.X
        logger_sctriangulate.info('running Scrublet may take several minutes')
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



    def compute_metrics(self,parallel=True):
        if parallel:
            cores1 = len(self.query)  # make sure to request same numeber of cores as the length of query list
            cores2 = mp.cpu_count()
            cores = min(cores1,cores2)
            logger_sctriangulate.info('Spawn to {} processes'.format(cores))
            pool = mp.Pool(processes=cores)
            self._to_sparse()
            raw_results = [pool.apply_async(each_key_run,args=(self,key)) for key in self.query]
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
            if self.reference not in self.query:
                self.cluster[self.reference] = self.adata.obs[self.reference].unique() 
            subprocess.run(['rm','-r','{}'.format(os.path.join(self.dir,'scTriangulate_local_mode_enrichr/'))])
            self._to_sparse()

        else:
            logger_sctriangulate.info('choosing to compute metrics sequentially')
            for key in self.query:
                collect = each_key_run(self,key)
                key = collect['key']
                for metric in self.metrics + list(self.add_metrics.keys()):
                    self.adata.obs['{}@{}'.format(metric,key)] = collect['col_{}'.format(metric)]
                self.score[key] = collect['score_info']
                self.cluster[key] = collect['cluster_info']  
                self._add_to_uns('confusion_reassign',key,collect)
                self._add_to_uns('confusion_sccaf',key,collect)
                self._add_to_uns('marker_genes',key,collect)
                self._add_to_uns('exclusive_genes',key,collect)
            if self.reference not in self.query:
                self.cluster[self.reference] = self.adata.obs[self.reference].unique() 
            subprocess.run(['rm','-r','{}'.format(os.path.join(self.dir,'scTriangulate_local_mode_enrichr/'))])
            self._to_sparse()
            

    def penalize_artifact(self,mode,stamps=None):
        '''void mode is to set stamp position to 0, stamp is like {leiden1:5}'''
        if mode == 'void':
            obs = self.adata.obs
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
            






    def compute_shapley(self,parallel=True):
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
            cores = mp.cpu_count()
            # split the obs and data, based on cell axis
            obs = self.adata.obs
            obs_index = np.arange(obs.shape[0])
            sub_indices = np.array_split(obs_index,cores)
            sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_obs, sub_obs, sub_obs]
            sub_datas = [data[:,sub_index,:] for sub_index in sub_indices]  # [sub_data,sub_data,....]
            pool = mp.Pool(processes=cores)
            logger_sctriangulate.info('spawn {} sub processes for shapley computing'.format(cores))
            raw_results = [pool.apply_async(func=run_shapley,args=(sub_obs[i],self.query,self.reference,self.size_dict,sub_datas[i])) for i in range(len(sub_obs))]
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
            cores = mp.cpu_count()
            sub_indices = np.array_split(obs_index,cores)  # indices for each chunk [(0,1,2...),(56,57,58...),(),....]
            sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_df,sub_df,...]
            pool = mp.Pool(processes=cores)
            logger_sctriangulate.info('spawn {} sub processes for getting raw sctriangulate result'.format(cores))
            r = pool.map_async(run_assign,sub_obs)
            pool.close()
            pool.join()
            results = r.get()  # [sub_obs,sub_obs...]
            obs = pd.concat(results)
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


    def pruning(self,method='reference',parallel=True):
        if parallel:
            if method == 'reference':
                obs = reference_pruning(self.adata.obs,self.reference,self.size_dict)
                self.adata.obs = obs

            elif method == 'reassign':
                obs, invalid = reassign_pruning(self)
                self.adata.obs = obs
                self.invalid = invalid
        self._prefixing(col='pruned')

        # finally, generate a celltype sheet
        obs = self.adata.obs
        with open(os.path.join(self.dir,'celltype.txt'),'w') as f:
            f.write('reference\tcell_cluster\tchoice\n')
            for ref,grouped_df in obs.groupby(by=self.reference):
                unique = grouped_df['pruned'].unique()
                for reassign in unique:
                    f.write('{}\t{}\n'.format(self.reference + '@' + ref,reassign))



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
        

    def plot_umap(self,col,kind='category',save=True):
        # col means which column in obs to draw umap on
        if kind == 'category':
            fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
            sc.pl.umap(self.adata,color=col,frameon=False,ax=ax[0])
            sc.pl.umap(self.adata,color=col,frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
            if save:
                plt.savefig(os.path.join(self.dir,'umap_sctriangulate_{}.png'.format(col)),bbox_inches='tight')
                plt.close()
        elif kind == 'continuous':
            sc.pl.umap(self.adata,color=col,frameon=False,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'umap_sctriangulate_{}.png'.format(col)),bbox_inches='tight')
                plt.close()

    def plot_confusion(self,name,key,save,**kwargs):
        df = self.uns[name][key]
        df = df.apply(func=lambda x:x/x.sum(),axis=1)
        sns.heatmap(df,cmap=scphere_cmap,**kwargs)  
        if save:
            plt.savefig(os.path.join(self.dir,'confusion_{}_{}.png'.format(name,key)),bbox_inches='tight')
            plt.close()
    
    def plot_cluster_feature(self,key,cluster,feature,enrichment_type='enrichr',save=False):
        if feature == 'enrichment':
            fig,ax = plt.subplots()
            a = self.uns['marker_genes'][key].loc[cluster,:][enrichment_type]
            ax.barh(y=np.arange(len(a)),width=[item for item in a.values()],color='#FF9A91')
            ax.set_yticks(np.arange(len(a)))
            ax.set_yticklabels([item for item in a.keys()])
            ax.set_title('Marker gene enrichment')
            ax.set_xlabel('-Log10(adjusted_pval)')
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_enrichment.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()
        elif feature == 'marker_gene':
            a = self.uns['marker_genes'][key].loc[cluster,:]['purify']
            top = a[:10]
            # change cmap a bit
            sc.pl.umap(self.adata,color=top,ncols=5,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_marker_umap.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()
        elif feature == 'exclusive_gene':
            a = self.uns['exclusive_genes'][key][cluster]  # self.uns['exclusive_genes'][key] is a pd.Series
            a = list(a.keys())
            top = a[:10]
            sc.pl.umap(self.adata,color=top,ncols=5,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_exclusive_umap.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()
        elif feature == 'location':
            col = [1 if item == str(cluster) else 0 for item in self.adata.obs[key]]
            self.adata.obs['tmp_plot'] = col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap('YlOrRd'),vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_location_umap.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()

    def plot_heterogeneity(self,cluster,style,save=False,genes=None): 
        adata_s = self.adata[self.adata.obs[self.reference]==cluster,:]  
        # remove prior color stamps
        tmp = adata_s.uns
        tmp.pop('prefixed_colors',None)
        adata_s.uns = tmp

        if style == 'build':  # draw umap and heatmap

            # umap
            fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0.5},figsize=(5,10))
            # ax1
            sc.pl.umap(adata_s,color=['prefixed'],ax=axes[0])
            # ax2
            col = [1 if item == str(cluster) else 0 for item in self.adata.obs[self.reference]]
            self.adata.obs['tmp_plot'] = col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap('YlOrRd'),vmin=1e-5,ax=axes[1])
            if save:
                plt.savefig(os.path.join(self.dir,'{}_heterogeneity_{}.png'.format(cluster,'umap')),bbox_inches='tight')
                plt.close()

            # heatmap
            tmp = adata_s.uns
            tmp.pop('rank_genes_groups',None)
            adata_s.uns = tmp

            if len(adata_s.obs['prefixed'].unique()) == 1: # it is already unique
                logger_sctriangulate.info('{0} entirely being assigned to one type, no need to do DE'.format(cluster))
                return None
            else:
                sc.tl.rank_genes_groups(adata_s,groupby='prefixed')
                adata_s = filter_DE_genes(adata_s,self.species,self.criterion)
                number_of_groups = len(adata_s.obs['prefixed'].unique())
                genes_to_pick = 50 // number_of_groups
                sc.pl.rank_genes_groups_heatmap(adata_s,n_genes=genes_to_pick,swap_axes=True,key='rank_genes_gruops_filtered')
                if save:
                    plt.savefig(os.path.join(self.dir,'{}_heterogeneity_{}.png'.format(cluster,'heatmap')),bbox_inches='tight')
                    plt.close()


        elif style == 'umap':
            fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0.5},figsize=(5,10))
            # ax1
            sc.pl.umap(adata_s,color=['prefixed'],ax=axes[0])
            # ax2
            col = [1 if item == str(cluster) else 0 for item in self.adata.obs[self.reference]]
            self.adata.obs['tmp_plot'] = col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=bg_greyed_cmap('YlOrRd'),vmin=1e-5,ax=axes[1])
            if save:
                plt.savefig(os.path.join(self.dir,'{}_heterogeneity_{}.png'.format(cluster,style)),bbox_inches='tight')
                plt.close()

        elif style == 'heatmap':
            tmp = adata_s.uns
            tmp.pop('rank_genes_groups',None)
            adata_s.uns = tmp
            if len(adata_s.obs['prefixed'].unique()) == 1: # it is already unique
                logger_sctriangulate.info('{0} entirely being assigned to one type, no need to do DE'.format(cluster))
                return None
            else:
                sc.tl.rank_genes_groups(adata_s,groupby='prefixed')
                adata_s = filter_DE_genes(adata_s,self.species,self.criterion)
                number_of_groups = len(adata_s.obs['prefixed'].unique())
                genes_to_pick = 50 // number_of_groups
                sc.pl.rank_genes_groups_heatmap(adata_s,n_genes=genes_to_pick,swap_axes=True,key='rank_genes_gruops_filtered')
                if save:
                    plt.savefig(os.path.join(self.dir,'{}_heterogeneity_{}.png'.format(cluster,style)),bbox_inches='tight')
                    plt.close()

        elif style == 'violin':
            sc.pl.violin(adata_s,genes,groupby='prefixed')
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}.png'.format(cluster,genes,style)),bbox_inches='tight')
                plt.close()
                
        elif style == 'cellxgene':
            adata_s.write(os.path.join(self.dir,'{}_heterogeneity_{}.h5ad'.format(cluster,style)))


        elif style == 'sankey':
            try:
                import plotly.graph_objects as go
                import kaleido
            except:
                logger_sctriangulate.warning('no plotly or kaleido library, fall back to matplotlib sankey plot')
                # processing the obs
                df = pd.DataFrame()
                df['ref'] = ['ref'+':'+item.split('|')[0] for item in adata_s.obs['prefixed']]   # ref:gs@ERP4
                df['query'] = [item.split('|')[1].split('@')[0] for item in adata_s.obs['prefixed']]  # leiden1
                df['cluster'] = [item.split('|')[1] for item in adata_s.obs['prefixed']]  # leiden1@5
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
                    plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}.png'.format(cluster,genes,style)),bbox_inches='tight')
                    plt.close()

            else:
                df = pd.DataFrame()
                df['ref'] = ['ref'+':'+item.split('|')[0] for item in adata_s.obs['prefixed']]   # ref:gs@ERP4
                df['query'] = [item.split('|')[1].split('@')[0] for item in adata_s.obs['prefixed']]  # leiden1
                df['cluster'] = [item.split('|')[1] for item in adata_s.obs['prefixed']]  # leiden1@5

                unique_ref = df['ref'].unique().tolist() # not lexicographically sorted, only one
                unique_query = df['query'].unique().tolist()  # not lexicographically sorted
                unique_cluster = df['cluster'].unique().tolist() # not lexicographically sorted

                # get node label and node color
                node_label = unique_ref + unique_query + unique_cluster
                from matplotlib import cm,colors
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
                link_color = [node_color[i] for i in link_source]

                print(node_label,node_color,link_source,link_target,link_value,link_color)

                # start to draw using plotly and save using kaleido
                node_plotly = dict(pad = 15, thickness = 20,line = dict(color = "black", width = 0.5),label = node_label,color = node_color)
                link_plotly = dict(source=link_source,target=link_target,value=link_value,color=link_color)
                fig = go.Figure(data=[go.Sankey(node = node_plotly,link = link_plotly)])
                fig.update_layout(title_text="{}_heterogeneity_{}".format(cluster,style), font_size=10)
                fig.write_image(os.path.join(self.dir,'{}_heterogeneity_{}.png'.format(cluster,style)))


    def plot_circular_barplot(self,col,save=False):
        # col can be 'raw' or 'pruned'
        obs = copy.deepcopy(self.adata.obs)
        reference = self.reference
        obs['value'] = np.full(shape=obs.shape[0], fill_value=1)
        obs = obs.loc[:, [reference, col, 'value']]
        obs4plot = obs.groupby(by=[reference, col])['value'].sum().reset_index()    
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
                    title='Reference:{}'.format(reference),borderaxespad=4)
        

        if save:
            plt.savefig(os.path.join(self.dir,'sctri_circular_barplot_{}.png'.format(col)),bbox_inches='tight')
            plt.close()

        



            

    def _atomic_viewer_figure(self,key):
        for cluster in self.cluster[key]:
            self.plot_cluster_feature(key,cluster,'enrichment','enrichr',True)
            self.plot_cluster_feature(key,cluster,'marker_gene','enrichr',True)
            self.plot_cluster_feature(key,cluster,'exclusive_gene','enrichr',True)
            self.plot_cluster_feature(key,cluster,'location','enrichr',True)

    def _atomic_viewer_hetero(self):
        for cluster in self.cluster[self.reference]:
            self.plot_heterogeneity(cluster,'build',True)



    def building_viewer_fig(self,parallel=True):
        logger_sctriangulate.info('Building viewer requires generating all the necessary figures, may take several minutes')
        # create a folder to store all the figures
        if not os.path.exists(os.path.join(self.dir,'figure4viewer')):
            os.mkdir(os.path.join(self.dir,'figure4viewer'))
        ori_dir = self.dir
        new_dir = os.path.join(self.dir,'figure4viewer')
        self.dir = new_dir
        # generate all the figures
        '''doublet plot'''
        self.plot_umap('doublet_scores','continuous',True)

        if platform.system() == 'Linux':    # can parallelize
            '''heterogeneity'''
            logger_sctriangulate.info('spawn 1 sub process for inspection figure genearation')
            p = mp.Process(target=self._atomic_viewer_hetero)
            p.start()
            '''other figures'''
            cores = mp.cpu_count()
            pool = mp.Pool(processes=cores)
            logger_sctriangulate.info('spawn {} sub processes for viewer figure generation'.format(cores))
            raw_results = [pool.apply_async(func=self._atomic_viewer_figure,args=(key,)) for key in self.cluster.keys()]
            p.join()
            pool.close()
            pool.join()
        else:                               # Windows and Darwin can not parallelize if plotting
            '''heterogeneity'''
            self._atomic_viewer_hetero()
            '''other figures'''
            for key in self.cluster.keys():
                self._atomic_viewer_figure(key)

        self.dir = ori_dir

    def building_viewer_html(self):
        with open(os.path.join(self.dir,'figure4viewer','viewer.html'),'w') as f:
            f.write(to_html(self.cluster,self.score,self.total_metrics))
        with open(os.path.join(self.dir,'figure4viewer','inspection.html'),'w') as f:
            f.write(inspection_html(self.cluster,self.reference))
        
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/viewer.js'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/viewer.css'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/inspection.js'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/inspection.css'),os.path.join(self.dir,'figure4viewer')))


        
        


# ancillary functions for main class
def penalize_artifact_void(obs,query,stamps,metrics):
    for stamp in stamps:
        metrics_cols = obs.loc[:,[item2+'@'+item1 for item1 in query for item2 in metrics]]
        cluster_cols = obs.loc[:,query]
        df = cluster_cols.apply(func=lambda x:pd.Series(data=[x.name+'@'+str(item) for item in x],name=x.name),axis=0)
        df_repeat = pd.DataFrame(np.repeat(df.values,len(metrics),axis=1))
        truth = pd.DataFrame(data=(df_repeat == stamp).values,index=metrics_cols.index,columns=metrics_cols.columns)
        tmp = metrics_cols.mask(truth,0)
        obs.loc[:,[item2+'@'+item1 for item1 in query for item2 in metrics]] = tmp
    return obs



def each_key_run(sctri,key):
    folder = sctri.dir
    adata = sctri.adata
    species = sctri.species
    criterion = sctri.criterion
    metrics = sctri.metrics
    add_metrics = sctri.add_metrics
    total_metrics = sctri.total_metrics

    try:
        assert issparse(adata.X) == False
    except AssertionError:
        adata.X = adata.X.toarray()  

    # remove cluster that only have 1 cell, for DE analysis
    adata_to_compute = check_filter_single_cluster(adata,key)  

    # a dynamically named dict
    cluster_to_metric = {}
    '''marker gene'''
    marker_genes = marker_gene(adata_to_compute,key,species,criterion,folder)
    logger_sctriangulate.info('Process {}, for {}, finished marker genes finding'.format(os.getpid(),key))
    '''reassign score'''
    cluster_to_metric['cluster_to_reassign'], confusion_reassign = reassign_score(adata_to_compute,key,marker_genes)
    logger_sctriangulate.info('Process {}, for {}, finished reassign score computing'.format(os.getpid(),key))
    '''tfidf10 score'''
    cluster_to_metric['cluster_to_tfidf10'], exclusive_genes = tf_idf10_for_cluster(adata_to_compute,key,species,criterion)
    logger_sctriangulate.info('Process {}, for {}, finished tfidf score computing'.format(os.getpid(),key))
    '''SCCAF score'''
    cluster_to_metric['cluster_to_SCCAF'], confusion_sccaf = SCCAF_score(adata_to_compute,key, species, criterion)
    logger_sctriangulate.info('Process {}, for {}, finished SCCAF score computing'.format(os.getpid(),key))
    '''doublet score'''
    cluster_to_metric['cluster_to_doublet'] = doublet_compute(adata_to_compute,key)
    logger_sctriangulate.info('Process {}, for {}, finished doublet score assigning'.format(os.getpid(),key))
    '''added other scores'''
    for metric,func in add_metrics.items():
        cluster_to_metric['cluster_to_{}'.format(metric)] = func(adata_to_compute,key,species,criterion)
        logging.info('Process {}, for {}, finished {} score computing'.format(os.getpid(),key,metric))


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

    return collect


def run_shapley(obs,query,reference,size_dict,data):
    logger_sctriangulate.info('process {} need to process {} cells for shapley computing'.format(os.getpid(),data.shape[1]))
    final = []
    intermediate = []
    for i in range(data.shape[1]):
        layer = data[:,i,:]
        result = []
        for j in range(layer.shape[0]):
            result.append(shapley_value(j,layer))
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

def filter_DE_genes(adata,species,criterion):
    de_gene = pd.DataFrame.from_records(adata.uns['rank_genes_groups']['names']) #column use field name, index is none by default, so incremental int value
    artifact = set(read_artifact_genes(species,criterion).index)
    de_gene.mask(de_gene.isin(artifact),inplace=True)
    adata.uns['rank_genes_gruops_filtered'] = adata.uns['rank_genes_groups'].copy()
    adata.uns['rank_genes_gruops_filtered']['names'] = de_gene.to_records(index=False)
    return adata












