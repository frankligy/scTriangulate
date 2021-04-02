#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from ast import literal_eval



def draw_enrich_plots(key):
    # for example, key=gs
    m = pd.read_csv('./scTriangulate_result/marker_{}.txt'.format(key),sep='\t',index_col=0)
    for cluster in m.index:
        a = literal_eval(m.loc[cluster,:]['enrichr'])
        fig,ax = plt.subplots()
        ax.barh(y=np.arange(len(a)),width=[item for item in a.values()],color='#FF9A91')
        ax.set_yticks(np.arange(len(a)))
        ax.set_yticklabels([item for item in a.keys()])
        ax.set_title('Marker gene enrichment')
        ax.set_xlabel('-Log10(adjusted_pval)')
        plt.savefig('./scTriangulate_diagnose/{0}_{1}_enrichment.png'.format(key,cluster),bbox_inches='tight')
        plt.close()

def draw_umap(adata,key):
    # for example, key=gs, cluster='Mono'
    m = pd.read_csv('./scTriangulate_result/marker_{}.txt'.format(key),sep='\t',index_col=0)
    for cluster in m.index:
        a = literal_eval(m.loc[cluster,:]['purify'])
        top = a[:10]
        sc.pl.umap(adata,color=top,ncols=5)
        plt.savefig('./scTriangulate_diagnose/{0}_{1}_marker_umap.png'.format(key,cluster),bbox_inches='tight')
        plt.close()

    e = pd.read_csv('./scTriangulate_result/exclusive_gene_{}.txt'.format(key),sep='\t',index_col=0)
    for cluster in e.index:
        a = literal_eval(e.loc[cluster,:]['genes'])
        top = a[:10]
        sc.pl.umap(adata,color=top,ncols=5)
        plt.savefig('./scTriangulate_diagnose/{0}_{1}_exclusive_umap.png'.format(key,cluster),bbox_inches='tight')
        plt.close()



