#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import anndata as ad
from scipy.stats import gmean
import pickle
from metrics import *

def check_filter_single_cluster(adata,key):
    vc = adata.obs[key].astype('str').value_counts()
    exclude_clusters= vc.loc[vc==1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    adata_valid = adata[truth,:]
    return adata_valid

adata = sc.read('scTriangulate_present/user_choice.h5ad')
adata.X = adata.X.toarray()

# all_reassign = []
# all_tfidf = []
# all_SCCAF = []
# for i in range(5):
#     random = np.random.choice(np.arange(15),size=(adata.X.shape[0],)).astype('str')
#     adata.obs['random'] = random
#     adata.obs['random'] = adata.obs['random'].astype('category')
#     key = 'random'
#     adata_to_compute = check_filter_single_cluster(adata,key)  # be a view
#     print('finished filtering')
#     result = marker_gene(adata_to_compute,key=key)
#     print('finished marker gene computing')
#     cluster_to_accuracy_random = reassign_score(adata_to_compute,key,result)
#     print('finished reassign score')
#     cluster_to_tfidf_random = tf_idf_for_cluster(adata_to_compute,key)
#     print('finished tfidf')
#     cluster_to_SCCAF_random = SCCAF_score(adata_to_compute,key)
#     print('finished SCCAF')

#     all_reassign.extend(list(cluster_to_accuracy_random.values()))
#     all_tfidf.extend(list(cluster_to_tfidf_random.values()))
#     all_SCCAF.extend(list(cluster_to_SCCAF_random.values()))

# data = pd.DataFrame({'all_reassign':all_reassign,'all_tfidf':all_tfidf,'all_SCCAF':all_SCCAF})
# sns.kdeplot(data=data,fill=True)
# plt.savefig('./sanity_check.png',bbox_inches='tight')
# plt.close()

# with open('random.p','wb') as f:
#     pickle.dump(data,f)

with open('random.p','rb') as f:
    data = pickle.load(f)

reassign_real = adata.obs['reassign@leiden1'].unique()
tfidf_real = np.concatenate([adata.obs['tfidf@leiden1'].unique(),adata.obs['tfidf@leiden2'].unique()])
SCCAF_real = adata.obs['SCCAF@leiden1'].unique()

sns.kdeplot(data=data,x='all_reassign',fill=True,color='blue')
sns.kdeplot(data=pd.DataFrame({'real_reassign':reassign_real}),fill=True,color='orange')
plt.savefig('./reassign.png',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
sns.kdeplot(data=data,x='all_tfidf',fill=True,color='blue',ax=ax)
sns.kdeplot(data=pd.DataFrame({'real_tfidf':tfidf_real}),fill=True,color='orange',ax=ax)
axins = ax.inset_axes([0.5,0.5,0.3,0.3])
axins.set_ylim([0,1.5])
sns.kdeplot(data=pd.DataFrame({'real_tfidf':tfidf_real}),fill=True,color='orange',ax=axins)
plt.savefig('./tfidf.png',bbox_inches='tight')
plt.close()

sns.kdeplot(data=data,x='all_SCCAF',fill=True,color='blue')
sns.kdeplot(data=pd.DataFrame({'real_SCCAF':SCCAF_real}),fill=True,color='orange')
plt.savefig('./SCCAF.png',bbox_inches='tight')
plt.close()
























