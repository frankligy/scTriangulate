
import pandas as pd
import numpy as np

result = pd.read_csv('/Volumes/salomonis2/LabFiles/Frank-Li/scTriangulate/pbmc10k/criterion2/scTriangulate_present/shapley_annotation.txt',
                     sep='\t',index_col=0)

from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, homogeneity_completeness_v_measure, v_measure_score
from sklearn.metrics import normalized_mutual_info_score
azimuth = LabelEncoder().fit_transform(result['azimuth'].values)
leiden1 = LabelEncoder().fit_transform(result['leiden1'].values)
leiden2 = LabelEncoder().fit_transform(result['leiden2'].values)
leiden3 = LabelEncoder().fit_transform(result['leiden3'].values)
final = LabelEncoder().fit_transform(result['reassign_prefix'].values)
user = LabelEncoder().fit_transform(result['choice'].values)

for test in [leiden1,leiden2,leiden3,final,azimuth]:
    print('ARI:{}'.format(adjusted_rand_score(azimuth,test)))
    print('AMI:{}'.format(adjusted_mutual_info_score(azimuth, test)))
    print('V:{}'.format(v_measure_score(azimuth, test)))
    print('NMI:{}'.format(normalized_mutual_info_score(azimuth,test)))