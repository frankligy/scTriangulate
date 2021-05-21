import sys
import os
import copy
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import seaborn as sns
from scipy.sparse import issparse,csr_matrix
import multiprocessing as mp
import platform
import logging
logging.basicConfig(format='%(asctime)s - %(message)s',level=logging.ERROR)

import scanpy as sc
import anndata as ad
import gseapy as gp
import scrublet as scr

from .shapley import *
from .metrics import *
from .viewer import *
from .prune import *



# define ScTriangulate Object
class ScTriangulate(object):

    def __init__(self,dir,adata,query,reference,species,criterion=2):
        self.dir = dir
        self.adata = adata
        self.query = query
        self.reference = reference
        self.species = species
        self.criterion = criterion
        self.score = {}
        self.cluster = {}
        self.uns = {}
        self._special_cmap()

    def __str__(self):  # when you print(instance) in REPL
        return 'ScTriangualate Object:\nWorking directory is {0}\nQuery Annotation is {1}\nReference Annotation is {2}\n'\
            'Species is {3}\nCriterion is {4}\nScores slot contains {5}\nCluster slot contains {6}\nUns slot contains {7}\n'\
            'cmap contains: {8}'.format(self.dir, self.query,self.reference,self.species,self.criterion,list(self.score.keys()),
            list(self.cluster.keys()),list(self.uns.keys()),list(self.cmap.keys()))

    def __repr__(self):  # when you type the instance in REPL
        return 'ScTriangualate Object:\nWorking directory is {0}\nQuery Annotation is {1}\nReference Annotation is {2}\n'\
            'Species is {3}\nCriterion is {4}\nScores slot contains {5}\nCluster slot contains {6}\nUns slot contains {7}\n'\
            'cmap contains: {8}'.format(self.dir, self.query,self.reference,self.species,self.criterion,list(self.score.keys()),
            list(self.cluster.keys()),list(self.uns.keys()),list(self.cmap.keys()))


    def _to_dense(self):
        self.adata.X = self.adata.X.toarray() 
        
    def _to_sparse(self):
        self.adata.X = csr_matrix(self.adata.X)

    def _special_cmap(self):
        self.cmap = {}
        cmap = copy.copy(cm.get_cmap('viridis'))
        cmap.set_under('lightgrey')
        self.cmap['viridis'] = cmap
        cmap = copy.copy(cm.get_cmap('YlOrRd'))
        cmap.set_under('lightgrey')
        self.cmap['YlOrRd'] = cmap


    def doublet_predict(self):
        if issparse(self.adata.X):
            self._to_dense()
        counts_matrix = self.adata.X
        logging.info('running Scrublet may take several minutes')
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
            logging.info('Spawn to {} processes'.format(cores))
            pool = mp.Pool(processes=cores)
            self._to_sparse()
            raw_results = [pool.apply_async(each_key_run_parallel,args=(self.adata,key,self.species,self.criterion)) for key in self.query]
            pool.close()
            pool.join()
            for collect in raw_results:
                collect = collect.get()
                key = collect['key']
                self.adata.obs['reassign@{}'.format(key)] = collect['col_reassign']
                self.adata.obs['tfidf@{}'.format(key)] = collect['col_tfidf']
                self.adata.obs['SCCAF@{}'.format(key)] = collect['col_SCCAF']
                self.score[key] = collect['score_info']
                self.cluster[key] = collect['cluster_info']  
                self._add_to_uns('confusion_reassign',key,collect)
                self._add_to_uns('confusion_sccaf',key,collect)
                self._add_to_uns('marker_genes',key,collect)
                self._add_to_uns('exclusive_genes',key,collect)
            if self.reference not in self.query:
                self.cluster[self.reference] = self.adata.obs[self.reference].unique() 
            os.system('rm -r "./scTriangulate_local_mode_enrichr/"')  # remove this temporary directory 

    def compute_shapley(self,parallel=True):
        if parallel:
            size_dict,size_list = get_size(self.adata.obs,self.query)
            self.size_dict = size_dict
            # compute shaley value
            score_colname = ['reassign','tfidf','SCCAF']
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
            sub_datas = np.array_split(data,cores,axis=1)  # [sub_data,sub_data,....]
            pool = mp.Pool(processes=cores)
            logging.info('spawn {} sub processes for shapley computing'.format(cores))
            raw_results = [pool.apply_async(func=run_shapley,args=(self.adata.obs,self.query,self.reference,self.size_dict,sub_data)) for sub_data in sub_datas]
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
            logging.info('spawn {} sub processes for getting raw sctriangulate result'.format(cores))
            r = pool.map_async(run_assign,sub_obs)
            pool.close()
            pool.join()
            results = r.get()  # [sub_obs,sub_obs...]
            obs = pd.concat(results)
            self.adata.obs = obs

    def pruning(self,parallel=True):
        if parallel:
            obs = reference_pruning(self.adata.obs,self.reference,self.size_dict)
            self.adata.obs = obs

            # also, prefix the pruned assignment with reference annotation
            col1 = self.adata.obs['pruned']
            col2 = self.adata.obs[self.reference]
            col = []
            for i in range(len(col1)):
                concat = self.reference + '@' + col2[i] + '|' + col1[i]
                col.append(concat)
            self.adata.obs['prefixed'] = col

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
                plt.savefig(os.path.join(self.dir,'umap_sctriangulate_{}.pdf'.format(col)),bbox_inches='tight')
                plt.close()
        elif kind == 'continuous':
            sc.pl.umap(self.adata,color=col,frameon=False,cmap=self.cmap['viridis'],vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'umap_sctriangulate_{}.pdf'.format(col)),bbox_inches='tight')
                plt.close()

    def plot_confusion(self,name,key,save,**kwargs):
        df = self.uns[name][key]
        df = df.apply(func=lambda x:x/x.sum(),axis=1)
        sns.heatmap(df,cmap='bwr',**kwargs)
        if save:
            plt.savefig(os.path.join(self.dir,'confusion_{}_{}.pdf'.format(name,key)),bbox_inches='tight')
            plt.close()
    
    def plot_cluster_feature(self,key,cluster,feature,save):
        if feature == 'enrichment':
            fig,ax = plt.subplots()
            a = self.uns['marker_genes'][key].loc[cluster,:]['enrichr']
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
            sc.pl.umap(self.adata,color=top,ncols=5,cmap=self.cmap['viridis'],vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_marker_umap.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()
        elif feature == 'exclusive_gene':
            a = self.uns['exclusive_genes'][key][cluster]  # self.uns['exclusive_genes'][key] is a pd.Series
            a = list(a.keys())
            top = a[:10]
            sc.pl.umap(self.adata,color=top,ncols=5,cmap=self.cmap['viridis'],vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_exclusive_umap.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()
        elif feature == 'location':
            col = [1 if item == str(cluster) else 0 for item in self.adata.obs[key]]
            self.adata.obs['tmp_plot'] = col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=self.cmap['YlOrRd'],vmin=1e-5)
            if save:
                plt.savefig(os.path.join(self.dir,'{0}_{1}_location_umap.png'.format(key,cluster)),bbox_inches='tight')
                plt.close()

    def plot_heterogeneity(self,cluster,style,save=False,genes=None):
        # cluster should be a valid cluster in self.reference
        adata_s = self.adata[self.adata.obs[self.reference]==cluster,:]
        if style == 'umap':
            fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0.5},figsize=(5,10))
            # ax1
            sc.pl.umap(adata_s,color=['prefixed'],ax=axes[0])
            # ax2
            col = [1 if item == str(cluster) else 0 for item in self.adata.obs[self.reference]]
            self.adata.obs['tmp_plot'] = col
            sc.pl.umap(self.adata,color='tmp_plot',cmap=self.cmap['YlOrRd'],vmin=1e-5,ax=axes[1])
            if save:
                plt.savefig(os.path.join(self.dir,'{}_heterogeneity_{}.png'.format(cluster,style)),bbox_inches='tight')
                plt.close()
        elif style == 'heatmap':
            if adata_s.uns.get('rank_genes_groups') != None:
                del adata_s.uns['rank_genes_groups']
            if len(adata_s.obs['prefixed'].unique()) == 1: # it is already unique
                logging.info('{0} entirely being assigned to one type, no need to do DE'.format(cluster))
                return None
            else:
                sc.tl.rank_genes_groups(adata_s,groupby='prefixed')
                adata_s = filter_DE_genes(adata_s,self.species,self.criterion)
                number_of_groups = len(adata_s.obs['prefixed'].unique())
                genes_to_pick = 50 // number_of_groups
                sc.pl.rank_genes_groups_heatmap(adata_s,n_genes=genes_to_pick,swap_axes=True,key='rank_genes_gruops_filtered')
                if save:
                    plt.savefig(os.path.join(self.dir,'{}_heterogeneity_{}.pdf'.format(cluster,style)),bbox_inches='tight')
                    plt.close()
        elif style == 'violin':
            sc.pl.violin(adata_s,genes,groupby='prefixed')
            if save:
                plt.savefig(os.path.join(self.dir,'{}_{}_heterogeneity_{}.png'.format(cluster,genes,style)),bbox_inches='tight')
                plt.close()
        elif style == 'cellxgene':
            adata_s.write(os.path.join(self.dir,'{}_heterogeneity_{}.h5ad'.format(cluster,style)))

    def _atomic_viewer_figure(self,key):
        for cluster in self.cluster[key]:
            self.plot_cluster_feature(key,cluster,'enrichment',True)
            self.plot_cluster_feature(key,cluster,'marker_gene',True)
            self.plot_cluster_feature(key,cluster,'exclusive_gene',True)
            self.plot_cluster_feature(key,cluster,'location',True)

    def _atomic_viewer_hetero(self):
        for cluster in self.cluster[self.reference]:
            self.plot_heterogeneity(cluster,'umap',True)
            self.plot_heterogeneity(cluster,'heatmap',True)



    def building_viewer_fig(self,parallel=True):
        logging.info('Building viewer requires generating all the necessary figures, may take several minutes')
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
            logging.info('spawn 1 sub process for inspection figure genearation')
            p = mp.Process(target=self._atomic_viewer_hetero)
            p.start()
            '''other figures'''
            cores = mp.cpu_count()
            pool = mp.Pool(processes=cores)
            logging.info('spawn {} sub processes for viewer figure generation'.format(cores))
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
            f.write(to_html(self.cluster,self.score))
        with open(os.path.join(self.dir,'figure4viewer','inspection.html'),'w') as f:
            f.write(inspection_html(self.cluster,self.reference))
        
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/viewer.js'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/viewer.css'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/inspection.js'),os.path.join(self.dir,'figure4viewer')))
        os.system('cp {} {}'.format(os.path.join(os.path.dirname(os.path.abspath(__file__)),'viewer/inspection.css'),os.path.join(self.dir,'figure4viewer')))


        
        


# ancillary functions for main class
def each_key_run_parallel(adata,key,species,criterion):
    adata.X = adata.X.toarray()
    assert issparse(adata.X) == False
    adata_to_compute = check_filter_single_cluster(adata,key)  
    marker_genes = marker_gene(adata_to_compute,key,species,criterion)
    logging.info('Process {}, for {}, finished marker genes finding'.format(os.getpid(),key))
    cluster_to_reassign, confusion_reassign = reassign_score(adata_to_compute,key,marker_genes)
    logging.info('Process {}, for {}, finished reassign score computing'.format(os.getpid(),key))
    cluster_to_tfidf, exclusive_genes = tf_idf_for_cluster(adata_to_compute,key,species,criterion)
    logging.info('Process {}, for {}, finished tfidf score computing'.format(os.getpid(),key))
    cluster_to_SCCAF, confusion_sccaf = SCCAF_score(adata_to_compute,key, species, criterion)
    logging.info('Process {}, for {}, finished SCCAF score computing'.format(os.getpid(),key))
    adata.X = csr_matrix(adata.X)
    assert issparse(adata.X) == True
    col_reassign = adata.obs[key].astype('str').map(cluster_to_reassign).fillna(0).values
    col_tfidf = adata.obs[key].astype('str').map(cluster_to_tfidf).fillna(0).values
    col_SCCAF = adata.obs[key].astype('str').map(cluster_to_SCCAF).fillna(0).values
    # add a doublet score to each cluster
    cluster_to_doublet = doublet_compute(adata_to_compute,key)
    score_info = [cluster_to_reassign,cluster_to_tfidf,cluster_to_SCCAF,cluster_to_doublet]
    cluster_info = list(cluster_to_reassign.keys())
    # all the intermediate results needed to be returned
    collect = {'key':key,
               'col_reassign':col_reassign,
               'col_tfidf':col_tfidf,
               'col_SCCAF':col_SCCAF,
               'score_info':score_info,
               'cluster_info':cluster_info,
               'marker_genes':marker_genes,
               'confusion_reassign':confusion_reassign,
               'exclusive_genes':exclusive_genes,
               'confusion_sccaf':confusion_sccaf,
               }
    return collect


def run_shapley(obs,query,reference,size_dict,data):
    logging.info('process {} need to process {} cells for shapley computing'.format(os.getpid(),data.shape[1]))
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
    logging.info('process {} need to process {} cells for raw sctriangulte result'.format(os.getpid(),obs.shape[0]))   
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

def plot_diagnose_figure(dir, adata, uns, cluster, cmap, key):
    for cluster in cluster[key]:
        # draw enrichment plot
        fig,ax = plt.subplots()
        a = uns['marker_genes'][key].loc[cluster,:]['enrichr']
        ax.barh(y=np.arange(len(a)),width=[item for item in a.values()],color='#FF9A91')
        ax.set_yticks(np.arange(len(a)))
        ax.set_yticklabels([item for item in a.keys()])
        ax.set_title('Marker gene enrichment')
        ax.set_xlabel('-Log10(adjusted_pval)')
        plt.savefig(os.path.join(dir,'{0}_{1}_enrichment.png'.format(key,cluster)),bbox_inches='tight')
        plt.close()

        # draw marker gene plot
        a = uns['marker_genes'][key].loc[cluster,:]['purify']
        top = a[:10]
        sc.pl.umap(adata,color=top,ncols=5,cmap=cmap['viridis'],vmin=1e-5)
        plt.savefig(os.path.join(dir,'{0}_{1}_marker_umap.png'.format(key,cluster)),bbox_inches='tight')
        plt.close()

        # draw exclusive gene plot
        a = uns['exclusive_genes'][key][cluster]  # self.uns['exclusive_genes'][key] is a pd.Series
        a = list(a.keys())
        top = a[:10]
        sc.pl.umap(adata,color=top,ncols=5,cmap=cmap['viridis'],vmin=1e-5)
        plt.savefig(os.path.join(dir,'{0}_{1}_exclusive_umap.png'.format(key,cluster)),bbox_inches='tight')
        plt.close()     

        # draw location plot
        col = [1 if item == str(cluster) else 0 for item in adata.obs[key]]
        adata.obs['tmp_plot'] = col
        sc.pl.umap(adata,color='tmp_plot',cmap=cmap['YlOrRd'],vmin=1e-5)
        plt.savefig(os.path.join(dir,'{0}_{1}_location_umap.png'.format(key,cluster)),bbox_inches='tight')
        plt.close()












