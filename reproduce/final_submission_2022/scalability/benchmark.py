#! /data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,sys
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# plot for effects on annotations

n_anno = [3,6,9,12]
time_run = [389,497,639,2090]
time_cpu = [1061,2511,6511,40515]

fig,ax = plt.subplots()
ax.plot(np.arange(len(n_anno)),[item/60 for item in time_run],linestyle='--',marker='o',label='run time')
ax.plot(np.arange(len(n_anno)),[item/60 for item in time_cpu],linestyle='--',marker='o',label='cpu time')
ax.set_xticks(np.arange(len(n_anno)))
ax.set_xticklabels(n_anno)
ax.set_xlabel('# annotations')
ax.set_ylabel('Running time (mins)')
ax.set_title('Effects of number of annotations on Running time -- CITE')
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('anno.pdf',bbox_inches='tight')
plt.close()


# plot for effects on cells

n_cells = [10000,20000,30000,40000]
time_run = [714,1341,1991,2456]

fig,ax = plt.subplots()
ax.plot(np.arange(len(n_cells)),[item/60 for item in time_run],linestyle='--',marker='o',label='run time')
ax.set_xticks(np.arange(len(n_cells)))
ax.set_xticklabels(n_cells)
ax.set_xlabel('# cells')
ax.set_ylabel('Running time (mins)')
ax.set_title('Effects of number of cells on Running time -- LUNG')
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('cells.pdf',bbox_inches='tight')
plt.close()


# plot for effects on 

n_features = [40000,60000,80000,100000]
time_run = [1233,1759,2547,3632]

fig,ax = plt.subplots()
ax.plot(np.arange(len(n_features)),[item/60 for item in time_run],linestyle='--',marker='o',label='run time')
ax.set_xticks(np.arange(len(n_features)))
ax.set_xticklabels(n_features)
ax.set_xlabel('# features')
ax.set_ylabel('Running time (mins)')
ax.set_title('Effects of number of features on Running time -- MULTIOME')
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
plt.savefig('features.pdf',bbox_inches='tight')
plt.close()















