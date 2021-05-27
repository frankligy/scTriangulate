#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import anndata as ad



def CLR_normalization(mat):
    '''matrix should be cell x ADT, expecting a ndarray'''
    from scipy.stats import gmean
    gmeans = gmean(mat+1,axis=1).reshape(-1,1)
    post = np.log(mat/gmeans + 1)
    return post

def total_normalization(mat,target=1e6):
    total = np.sum(mat,axis=1).reshape(-1,1)
    sf = total/target
    post = np.log(mat/sf + 1)
    return post

def GMM_normalization(mat):
    mat = total_normalization(mat)
    from sklearn.mixture import GaussianMixture
    model = GaussianMixture(n_components=2,random_state=0)
    model.fit(mat)
    means = model.means_  # (n_components,n_features)
    bg_index = np.argmin(means.mean(axis=1))
    bg_mean = means[bg_index,:].reshape(1,-1)
    post = mat - bg_mean
    return post

cite = pd.read_csv('./exp.Original-60-ADT-cKit.txt',sep='\t',index_col=0).drop(index='unmapped')
cite.rename(mapper=lambda x:x.split('-')[0],axis=0,inplace=True)
cite.rename(mapper=lambda x:x.replace('/','_'),axis=0,inplace=True)
cite.rename(mapper=lambda x:x.split('.')[0],axis=1,inplace=True)

cite.to_csv('ori_cite.txt',sep='\t')
clr_transform = pd.DataFrame(data=CLR_normalization(cite.values.T).T,index=cite.index,columns=cite.columns)
total_transform = pd.DataFrame(data=total_normalization(cite.values.T).T,index=cite.index,columns=cite.columns)
gmm_transform = pd.DataFrame(data=GMM_normalization(cite.values.T).T,index=cite.index,columns=cite.columns)

# clr_transform.to_csv('clr_cite.txt',sep='\t')
# total_transform.to_csv('total_cite.txt',sep='\t')
# gmm_transform.to_csv('gmm_cite.txt',sep='\t')

# plotting
for i,adt in enumerate(cite.index):
    fig,ax = plt.subplots()
    data = pd.DataFrame({'ori':cite.loc[adt,:].values,'clr':clr_transform.loc[adt,:].values,
                        'total':total_transform.loc[adt,:],'gmm':gmm_transform.loc[adt,:].values})
    sns.histplot(data=data,kde=True,ax=ax)
    ax.set_title('{}'.format(adt))
    ax.set_xlim([-10,20])
    plt.savefig('./dist/dist_{}.pdf'.format(adt),bbox_inches='tight')
    plt.close()

sys.exit('stop')







fig,axes = plt.subplots(nrows=16,ncols=5,figsize=(30,30),gridspec_kw={'wspace':0.3,'hspace':0.3})
axes = axes.flatten()
for i,adt in enumerate(cite.index):
    data = pd.DataFrame({'ori':cite.loc[adt,:].values,'clr':clr_transform.loc[adt,:].values,
                        'total':total_transform.loc[adt,:],'gmm':gmm_transform.loc[adt,:].values})
    sns.histplot(data=data,kde=True,ax=axes[i])
    axes[i].set_tile('{}'.format(adt))
plt.savefig('distribution.pdf',bbox_inches='tight')
plt.close()


     

    









