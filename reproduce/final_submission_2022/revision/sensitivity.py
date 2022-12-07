#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
import os,sys
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


metric = 'vmeasure'

data_dict = {
    'pbmc10k_scRNA':'/data/salomonis2/LabFiles/Frank-Li/scTriangulate/pbmc10k_qc/output_two',
    'TNC1_CITE':'/data/salomonis2/LabFiles/Frank-Li/scTriangulate/TNC1_qc/output_one',
    'TNC2_CITE':'/data/salomonis2/LabFiles/Frank-Li/scTriangulate/TNC2_qc/output_one',
    'Multiome':'/data/salomonis2/LabFiles/Frank-Li/scTriangulate/multiome_pbmc_qc/output_two',
    'TEA-Seq':'/data/salomonis2/LabFiles/Frank-Li/scTriangulate/TEA_qc/output_one'
}

# plot shapley_mode
options = ['pruned_rank_all_or_none','pruned_rank','pruned_shapley','pruned']
options_display = ['Rank(All_or_None)','Rank(classic)','Shapley(classic)','Shapley(All_or_None)-Default']
df_data = []
for ds,path in data_dict.items():
    full = os.path.join(path,'cluster_performance_shapley_mode.txt')
    tmp = pd.read_csv(full,sep='\t',index_col=0)
    df_data.append(tmp.loc[metric,options].values)
df = pd.DataFrame.from_records(df_data,index=list(data_dict.keys()),columns=options_display)
ax = df.plot.bar()
ax.get_legend().set_visible(False)
ax.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('sensitivity_shapley_mode_{}.pdf'.format(metric),bbox_inches='tight')
plt.close()


# plot bonus
options = ['pruned_bonus_0','pruned_bonus_0.005','pruned_bonus_0.01','pruned_bonus_0.05','pruned_bonus_0.1']
options_display = ['0','0.005','0.01','0.05','0.1']
df_data = []
for ds,path in data_dict.items():
    full = os.path.join(path,'cluster_performance_bonus.txt')
    tmp = pd.read_csv(full,sep='\t',index_col=0)
    df_data.append(tmp.loc[metric,options].values)
df = pd.DataFrame.from_records(df_data,index=list(data_dict.keys()),columns=options_display)
ax = df.plot.bar()
ax.get_legend().set_visible(False)
ax.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('sensitivity_bonus_{}.pdf'.format(metric),bbox_inches='tight')
plt.close()

# plot wfc
options = ['pruned_wfc_0','pruned_wfc_0.05','pruned_wfc_0.1','pruned_wfc_0.15','pruned_wfc_0.2','pruned_wfc_0.25','pruned_wfc_0.3','pruned_wfc_0.35','pruned_wfc_0.4','pruned_wfc_0.45','pruned_wfc_0.5','pruned_wfc_0.55','pruned_wfc_0.6']
options_display = ['0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']
df_data = []
for ds,path in data_dict.items():
    full = os.path.join(path,'cluster_performance_wfc.txt')
    tmp = pd.read_csv(full,sep='\t',index_col=0)
    df_data.append(tmp.loc[metric,options].values)
df = pd.DataFrame.from_records(df_data,index=list(data_dict.keys()),columns=options_display)
ax = df.plot.bar()
ax.get_legend().set_visible(False)
ax.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('sensitivity_wfc_{}.pdf'.format(metric),bbox_inches='tight')
plt.close()

# plot nabs
options = ['pruned_nabs_1','pruned_nabs_5','pruned_nabs_10','pruned_nabs_15','pruned_nabs_20','pruned_nabs_25','pruned_nabs_30','pruned_nabs_35','pruned_nabs_40','pruned_nabs_45','pruned_nabs_50']
options_display = ['1','5','10','15','20','25','30','35','40','45','50']
df_data = []
for ds,path in data_dict.items():
    full = os.path.join(path,'cluster_performance_nabs.txt')
    tmp = pd.read_csv(full,sep='\t',index_col=0)
    df_data.append(tmp.loc[metric,options].values)
df = pd.DataFrame.from_records(df_data,index=list(data_dict.keys()),columns=options_display)
ax = df.plot.bar()
ax.get_legend().set_visible(False)
ax.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('sensitivity_nabs_{}.pdf'.format(metric),bbox_inches='tight')
plt.close()

