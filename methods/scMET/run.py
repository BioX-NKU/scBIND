import pandas as pd
import scanpy as sc
import mofapy2
from mofapy2.run.entry_point import entry_point
import numpy as np
from sklearn.cluster import KMeans
import anndata as ad
import mofax as mfx
import os
import sys

h5admet = sys.argv[1]
h5adtot = sys.argv[2]
groundtruth = sys.argv[3]
n_clusters = int(sys.argv[4])
predcluster = sys.argv[5]
savepath = sys.argv[6]

hvf = pd.read_csv('hvf.csv', header=None, index_col=None)
hvf = hvf.iloc[:, 0].tolist()
adata_met = sc.read(h5admet)
adata_tot = sc.read(h5adtot)
adata = ad.AnnData(X=adata_met.X / adata_tot.X, obs=adata_met.obs.copy(), var = adata_met.var.copy())
adata = adata[:, hvf].copy()

data = [[adata.X]]
my_model = entry_point()
my_model.set_data_matrix(data, views_names=["expression"], samples_names=[adata.obs_names.tolist()], groups_names=["group0"])
my_model.set_model_options(factors=15)
my_model.set_train_options(seed=42)
my_model.build()
my_model.run()
my_model.save(outfile="mofa_model.hdf5")
model = mfx.mofa_model("mofa_model.hdf5")
cell_factors = model.get_factors()
model.close()

kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init='auto')
cluster_labels = kmeans.fit_predict(cell_factors)

df_clustering = pd.DataFrame()
df_clustering['scMET'] = cluster_labels
df_clustering.to_csv(predcluster, header=True, index=False)
true_label = pd.read_csv(groundtruth, sep='\t', header=0, index_col=None, compression='gzip').iloc[:, 1].values
print(cluster_labels)
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

cluster_metrics(true_label, cluster_labels)


