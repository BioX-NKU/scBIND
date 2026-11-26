from Step_2 import find_kopcluster
from Step_2 import find_kspcluster
from Step_2 import scM
from Step_2 import hc_pre
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import sys
import os

h5adpath = sys.argv[1]
groundtruth = sys.argv[2]
predcluster = sys.argv[3]
savepath = sys.argv[4]
n_clusters = sys.argv[5]

Cosine_Mat = pd.read_csv("./distance/dism_cosine.csv", header=None, index_col=None)
Hamming_Mat = pd.read_csv("./distance/dism_hamming.csv", header=None, index_col=None)
Pearson_Mat = pd.read_csv("./distance/dism_pearson.csv", header=None, index_col=None)
true_label = pd.read_csv(groundtruth, sep='\t', header=0, index_col=None, compression='gzip').iloc[:, 1].values

res_mat = scM(Cosine_Mat,Hamming_Mat,Pearson_Mat,n_clusters) ### the reconstructed similarity matrix 
C = hc_pre(res_mat,n_clusters) ### the resulting cell partitions

print(C)
adata = sc.read(h5adpath)
sorted_cell_names = sorted(adata.obs_names)
cell2major = dict(zip(adata.obs.index, C))
major_type_list = [cell2major[cell] for cell in sorted_cell_names]
major_type_arr = np.array(major_type_list)
df_clustering = pd.DataFrame()
df_clustering['scMelodyRegion-C'] = major_type_arr
df_clustering.to_csv(predcluster, header=True, index=False)
print(true_label)

from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import completeness_score

def cluster_metrics(target, pred):
    target = np.array(target)
    pred = np.array(pred)

    ari = adjusted_rand_score(target, pred)
    ami = adjusted_mutual_info_score(target, pred)
    nmi = normalized_mutual_info_score(target, pred)
    fmi = fowlkes_mallows_score(target, pred)
    comp = completeness_score(target, pred)
    homo = homogeneity_score(target, pred)
    print("")
    print(f"ari=%.6f" % (ari))
    print(f"ami=%.6f" % (ami))
    print(f"nmi=%.6f" % (nmi))
    print(f"fmi=%.6f" % (fmi))
    print(f"comp=%.6f" % (comp))
    print(f"homo=%.6f" % (homo))
    df = pd.DataFrame({'index': ['ari', 'ami', 'nmi', 'fmi', 'comp', 'homo'], 'value': [ari, ami, nmi, fmi, comp, homo]})
    df.to_csv(savepath, header=True, index=False, sep='\t')

cluster_metrics(true_label, major_type_arr)
