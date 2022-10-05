import pandas as pd
import squidpy as sq
import scanpy as sc
import os,sys
from sctriangulate.colors import *
from sctriangulate.preprocessing import *
from sctriangulate import *
from sctriangulate.spatial import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import f1_score, confusion_matrix, ConfusionMatrixDisplay
from tqdm import tqdm
from functools import reduce
from scipy.sparse import csr_matrix

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

os.chdir('/Users/ligk2e/Desktop/sitc')


def prepare_inputs(adata_sc,adata_sc_processed,adata_spatial,outdir='inputs',name=None,column=None):
    # sc, count mtx
    adata_to_mtx(adata_sc, outdir=os.path.join(outdir,'{}_sc_count_data'.format(name)))
    # sc, metadata
    adata_sc.obs['sample'] = np.full(adata_sc.shape[0], name)
    tmp = adata_sc.obs.loc[:, [column, 'sample']]
    tmp.index.name = 'index'
    tmp.to_csv(os.path.join(outdir,'{}_sc_meta.txt'.format(name)), sep='\t')
    tmp = adata_sc.obs.loc[:, [column]]
    tmp.index.name = 'index'
    tmp.to_csv(os.path.join(outdir,'{}_sc_meta_just_cell_type.txt'.format(name)), sep='\t')
    # spatial count and mtx and location
    adata_to_mtx(adata_spatial,outdir=os.path.join(outdir,'{}_spatial_count_data'.format(name)))
    tmp = adata_spatial.to_df().T
    tmp.index.name = 'index'
    tmp.to_csv(os.path.join(outdir,'{}_spatial_count.txt'.format(name)),sep='\t')
    tmp = pd.DataFrame(data=adata_spatial.obsm['spatial'], index=adata_spatial.obs_names, columns=['x', 'y'])
    tmp.index.name = 'index'
    tmp.to_csv(os.path.join(outdir,'{}_spatial_location.txt'.format(name)), sep='\t')
    # cell2location input
    cell2loc_input_adata_sc = adata_sc.copy()
    make_sure_adata_writable(cell2loc_input_adata_sc)
    cell2loc_input_adata_sc.X = csr_matrix(cell2loc_input_adata_sc.X)
    cell2loc_input_adata_sc.write(os.path.join(outdir,'cell2loc_input_adata_sc.h5ad'))
    cell2loc_input_adata_spatial = adata_spatial.copy()
    make_sure_adata_writable(cell2loc_input_adata_spatial)
    cell2loc_input_adata_spatial.X = csr_matrix(cell2loc_input_adata_spatial.X)
    cell2loc_input_adata_spatial.write(os.path.join(outdir,'cell2loc_input_adata_spatial.h5ad'))
    # seurat input
    seurat_spatial = adata_spatial.copy()
    seurat_spatial.X = csr_matrix(seurat_spatial.X)
    seurat_spatial.uns.pop('spatial')
    del seurat_spatial.obsm['spatial']
    seurat_spatial = make_sure_adata_writable(seurat_spatial)
    seurat_spatial.write(os.path.join(outdir,'seurat_spatial.h5ad'))
    seurat_sc = adata_sc.copy()
    seurat_sc.X = csr_matrix(seurat_sc.X)
    seurat_sc.write(os.path.join(outdir,'seurat_sc.h5ad'))
    # signature
    markers = pd.DataFrame.from_records(adata_sc_processed.uns['rank_genes_groups_filtered']['names'])
    n = 20
    marker_list = []
    for col in markers.columns:
        for i, gene in enumerate(markers[col]):
            if not pd.isnull(gene) and i < n:
                marker_list.append(gene)
    df_tmp = adata_sc_processed.to_df()
    df_tmp['celltype'] = adata_sc_processed.obs[column]
    sub_df_list = []
    ct_list = []
    for ct, sub_df in df_tmp.groupby('celltype'):
        sub_df_copy = sub_df.copy()
        sub_df_copy.drop(columns=['celltype'])
        sub_df_list.append(sub_df_copy.loc[:, marker_list].mean(axis=0))
        ct_list.append(ct)
    signatures = pd.concat(sub_df_list, axis=1)
    signatures.columns = ct_list
    signatures = signatures.loc[np.logical_not(signatures.index.duplicated()), :]
    signatures.to_csv(os.path.join(outdir,'{}_signatures.txt'.format(name)), sep='\t')


def run_ingest(adata_sc,adata_spatial,column,outdir):
    var_names = adata_sc.var_names.intersection(adata_spatial.var_names)
    adata_ref = adata_sc[:, var_names].copy()
    adata = adata_spatial[:, var_names].copy()
    def preprocess(adata):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=3000)
        adata = adata[:, adata.var['highly_variable']]
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        return adata
    adata_ref = preprocess(adata_ref)
    adata = preprocess(adata)
    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names].copy()
    adata = adata[:, var_names].copy()
    sc.tl.ingest(adata, adata_ref, obs=column)
    adata.obs['value'] = np.full(adata.shape[0],1)
    tmp = adata.obs.copy().reset_index()
    ingest_prop = pd.crosstab(index=tmp['index'],columns=tmp[column],values=tmp['value'],aggfunc=np.mean).fillna(value=0)
    missing_ct = list(set(adata_sc.obs[column].cat.categories).difference(ingest_prop.columns))
    missing_prop = pd.DataFrame(data=np.full(shape=(ingest_prop.shape[0],len(missing_ct)),fill_value=0),index=ingest_prop.index,columns=missing_ct)
    ingest_prop = pd.concat([ingest_prop,missing_prop],axis=1)
    ingest_prop.to_csv(os.path.join(outdir,'ingest_prop.txt'),sep='\t')

def run_scanorama(adata_sc,adata_spatial,column,outdir):
    import scanorama
    adata_sc.X = csr_matrix(adata_sc.X)
    adata_spatial.X = csr_matrix(adata_spatial.X)
    adata_sc.obs.index.name = None
    adata_spatial.obs.index.name = None
    adata_sc.var.index.name = None
    adata_spatial.var.index.name = None
    adatas = [adata_sc,adata_spatial]
    corrected = scanorama.correct_scanpy(adatas, return_dimred=True)
    adata_combined = ad.concat(corrected,axis=0,join='outer',merge='first')
    adata_combined_have_label = adata_combined[adata_combined.obs[column].notna(),:]
    adata_combined_no_label = adata_combined[adata_combined.obs[column].isna(),:]

    knn_X_train = adata_combined_have_label.obsm['X_scanorama']
    knn_X_test = adata_combined_no_label.obsm['X_scanorama']
    knn_Y_train = adata_combined_have_label.obs[column].values
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder().fit(knn_Y_train)
    knn_Y_train_encoded = le.transform(knn_Y_train)
    model = KNeighborsClassifier(n_neighbors=15,weights='distance')
    model.fit(knn_X_train,knn_Y_train_encoded)
    knn_Y_test = le.inverse_transform(model.predict(knn_X_test))
    adata_combined_no_label.obs[column] = knn_Y_test
    adata_combined_no_label.uns['spatial'] = adata_spatial.uns['spatial']
    adata_combined_no_label.obs['value'] = np.full(adata_combined_no_label.shape[0],1)
    tmp = adata_combined_no_label.obs.copy().reset_index()
    scanorama_prop = pd.crosstab(index=tmp['index'],columns=tmp[column],values=tmp['value'],aggfunc=np.mean).fillna(value=0)
    missing_ct = list(set(adata_sc.obs[column].cat.categories).difference(scanorama_prop.columns))
    missing_prop = pd.DataFrame(data=np.full(shape=(scanorama_prop.shape[0],len(missing_ct)),fill_value=0),index=scanorama_prop.index,columns=missing_ct)
    scanorama_prop = pd.concat([scanorama_prop,missing_prop],axis=1)
    scanorama_prop.to_csv(os.path.join(outdir,'scanoroma_prop.txt'),sep='\t')


def visualize_decon_result(decon_path,adata_spatial,sep='\t',output=None,markers_to_plot=None,outdir=None):
    decon = pd.read_csv(decon_path, sep=sep, index_col=0)
    common = decon.index.intersection(adata_spatial.obs_names)
    adata_vis = adata_spatial[common, :].copy()
    adata_vis.obsm['decon'] = decon.loc[common, :]
    sc.pp.normalize_total(adata_vis, target_sum=1e4)
    sc.pp.log1p(adata_vis)
    library_id = list(adata_spatial.uns['spatial'].keys())[0]
    types = ['B cells','Cancer Epithelial','Macrophage','Endothelial','T cells','fibroblast','Normal Epithelial','NK cells','DCs']
    sc.pl.spatial(sq.pl.extract(adata_vis, 'decon'),color=types+markers_to_plot, img_key='hires', scale_factor=adata_spatial.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef'])
    plt.savefig(os.path.join(outdir,'{}_vis_on_image.pdf'.format(output)), bbox_inches='tight')


def evaluate_decon_result_single(decon_path,adata_spatial,sep='\t',output=None,cutoff=0.05,outdir=None,plot_info=None):
    decon = pd.read_csv(decon_path, sep=sep, index_col=0)
    common = decon.index.intersection(adata_spatial.obs_names)
    adata_vis = adata_spatial[common, :].copy()
    adata_vis.obsm['decon'] = decon.loc[common, :]
    sc.pp.normalize_total(adata_vis, target_sum=1e4)
    sc.pp.log1p(adata_vis)
    # see plot_dict is ('Tcell','CD3')
    fig,ax = plt.subplots()
    cell,marker = plot_info
    value_cell = decon[cell].values
    value_marker = adata_vis[:,marker].X.toarray().squeeze()
    corr = pearsonr(value_cell, value_marker)[0]
    cutoff_cell = decon[cell].quantile(cutoff)
    cutoff_marker= np.quantile(adata_vis[:,marker].X.toarray().squeeze(),q=cutoff)
    sns.regplot(x=value_cell,y=value_marker,ax=ax)
    ax.axvline(cutoff_cell,linestyle='--',color='k')
    ax.axhline(cutoff_marker,linestyle='--',color='k')
    ax.set_xlabel('{} proportion'.format(cell))
    ax.set_ylabel('{}'.format(marker))
    ax.text(x=0.8,y=0.8,s='r={}'.format(round(corr,4)),transform=ax.transAxes)
    plt.savefig('single_{}_{}.pdf'.format(cell,marker),bbox_inches='tight')
    plt.close()



def evaluate_decon_result_merge(decon_path,adata_spatial,sep='\t',output=None,cutoff=0.05,outdir=None,markers_dict=None):
    decon = pd.read_csv(decon_path, sep=sep, index_col=0)
    common = decon.index.intersection(adata_spatial.obs_names)
    adata_vis = adata_spatial[common, :].copy()
    adata_vis.obsm['decon'] = decon.loc[common, :]
    sc.pp.normalize_total(adata_vis, target_sum=1e4)
    sc.pp.log1p(adata_vis)
    markers = reduce(lambda a,b:a+b,markers_dict.values())
    markers = list(set(markers))
    df_markers = adata_vis[:, markers].to_df()
    df_props = sq.pl.extract(adata_vis, 'decon').obs.loc[:, decon.columns]
    df = pd.concat([df_markers, df_props], axis=1)
    ncols = 4
    nrows = len(markers_dict) // ncols + 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.8})
    axes = axes.flatten()
    i = -1
    for k,v in markers_dict.items():
        i += 1
        df['average'] = df.loc[:,v].mean(axis=1).values
        sns.regplot(x=k,y='average',data=df,ax=axes[i],scatter_kws={'s':5})
        corr = pearsonr(df[k].values, df['average'].values)[0]
        cutoff_exp = df['average'].quantile(cutoff)
        cutoff_prop = df[k].quantile(cutoff)
        print('cutoff for {}: exp is {}, prop is {}'.format(k,cutoff_exp,cutoff_prop))
        y_true = np.where(df['average'].values > cutoff_exp, 1, 0)
        y_pred = np.where(df[k].values > cutoff_prop, 1, 0)
        try:
            f1 = f1_score(y_true, y_pred)
        except ValueError:  # all zero in y_true or y_pred
            f1 = 0
        axes[i].text(x=0.8,y=0.8,s='r={}\nf1={}'.format(round(corr,4),round(f1,4)),transform=axes[i].transAxes)
        axes[i].axvline(cutoff_prop,linestyle='--',color='k')
        axes[i].axhline(cutoff_exp,linestyle='--',color='k')
    plt.savefig(os.path.join(outdir,'{}_evaluate_merge.pdf'.format(output)), bbox_inches='tight')
    plt.close()

def build_the_df_merge(method,decon,adata_spatial,cutoff=0.05,markers_dict=None):
    decon = pd.read_csv(decon, sep='\t', index_col=0)
    common = decon.index.intersection(adata_spatial.obs_names)
    adata_vis = adata_spatial[common, :].copy()
    adata_vis.obsm['decon'] = decon.loc[common, :]
    sc.pp.normalize_total(adata_vis, target_sum=1e4)
    sc.pp.log1p(adata_vis)
    markers = reduce(lambda a,b:a+b,markers_dict.values())
    markers = list(set(markers))
    df_markers = adata_vis[:,markers].to_df()
    df_props = sq.pl.extract(adata_vis, 'decon').obs.loc[:,decon.columns]
    df = pd.concat([df_markers,df_props],axis=1)
    mapping = markers_dict
    data = []
    for k,v in mapping.items():
        df['average'] = df.loc[:,v].mean(axis=1).values
        corr_pearson = pearsonr(df[k].values, df['average'].values)[0]
        corr_spearman = spearmanr(df[k].values, df['average'].values)[0]
        cutoff_exp = df['average'].quantile(cutoff)
        cutoff_prop = df[k].quantile(cutoff)
        print('cutoff for {}: exp is {}, prop is {}'.format(k,cutoff_exp,cutoff_prop))
        y_true = np.where(df['average'].values > cutoff_exp, 1, 0)
        y_pred = np.where(df[k].values > cutoff_prop, 1, 0)
        try:
            f1 = f1_score(y_true, y_pred)
        except ValueError:  # all zero in y_true or y_pred
            f1 = 0
        data.append((method,'{}'.format(k),corr_pearson,corr_spearman,f1))
    tmp = pd.DataFrame.from_records(data,columns=['method','cell_type','pearsonr','spearmanr','f1_score'])
    return tmp


















