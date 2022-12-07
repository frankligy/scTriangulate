#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/ClusterEnsembles_env/bin/python3.7

import numpy as np
import ClusterEnsembles as CE
import pandas as pd
import sys

# read
obs = pd.read_csv('output_opt1_RN_two_new/to_bm_ensemble.txt',sep='\t',index_col=0)
query = ['Sun','Columbia','Kaminsky','Krasnow']
obs = obs.loc[:,query]

# encode
from sklearn.preprocessing import LabelEncoder
result = obs
labels = []
for r in obs.columns:
    encoded = LabelEncoder().fit_transform(result[r].values)
    labels.append(encoded)
labels = np.vstack(labels)

# Ensemble
#label_cspa = CE.cluster_ensembles(labels,solver='cspa')
label_hgpa = CE.cluster_ensembles(labels,solver='hgpa')
print('finished hgpa')
label_mcla = CE.cluster_ensembles(labels,solver='mcla')
print('finished mcla')
label_hbgf = CE.cluster_ensembles(labels,solver='hbgf')
print('finished hbgf')
#label_nmf = CE.cluster_ensembles(labels,solver='nmf')


# write
result = pd.DataFrame(data={'hgpa':label_hgpa,'mcla':label_mcla,'hbgf':label_hbgf,},index=obs.index)
result.to_csv('output_opt1_RN_two_new/from_bm_ensemble.txt',sep='\t')


