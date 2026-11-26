import episcanpy as epi
import numpy as np
import scanpy as sc
import sklearn
import pandas as pd
import sys
import os

predcluster = sys.argv[1]
savepath = sys.argv[2]

adata = sc.read_h5ad('constructed.h5ad')

sc.pp.neighbors(adata, n_neighbors=5, use_rep="X", metric="euclidean")
sc.tl.leiden(adata, resolution=1, key_added="leiden")

cluster_key='cell_type'
label_key='leiden'
print('Leiden')
print(adata.obs[cluster_key])
print(adata.obs[label_key])
df_clustering = pd.DataFrame()
df_clustering['MethSCAn'] = adata.obs[label_key]
df_clustering.to_csv(predcluster, header=True, index=False)
ARI = epi.tl.ARI(adata, cluster_key, label_key)
AMI = epi.tl.AMI(adata, cluster_key, label_key)
NMI = sklearn.metrics.normalized_mutual_info_score(adata.obs[cluster_key], adata.obs[label_key])
FMI = sklearn.metrics.fowlkes_mallows_score(adata.obs[cluster_key], adata.obs[label_key])
Comp = sklearn.metrics.completeness_score(adata.obs[cluster_key], adata.obs[label_key])
Homo = epi.tl.homogeneity(adata, cluster_key, label_key)
print('ARI:%.6f\tAMI:%.6f\tNMI:%.6f\tFMI:%.6f\tComp:%.6f\tHomoL%.6f'%(ARI, AMI, NMI, FMI, Comp, Homo))
df = pd.DataFrame({'index': ['ari', 'ami', 'nmi', 'fmi', 'comp', 'homo'], 'value': [ARI, AMI, NMI, FMI, Comp, Homo]})
df.to_csv(savepath, header=True, index=False, sep='\t')

