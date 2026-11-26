import snapatac2 as snap
import scanpy as sc
import anndata as ad
import numpy as np
from scipy import sparse
import pandas as pd
import sys

groundtruth = sys.argv[1]
h5adpath = sys.argv[2]
predcluster = sys.argv[3]
savepath = sys.argv[4]
n_clusters = int(sys.argv[5])

true_label = pd.read_csv(groundtruth, sep='\t', header=0, index_col=None, compression='gzip').iloc[:, 1].values
adata = sc.read(h5adpath)
adata = adata[:, ~np.all(np.isnan(adata.X), axis=0)].copy()
# Step 2: Calculate the median for each column, ignoring NaNs
medians = np.nanmedian(adata.X, axis=0)
# Step 3: Replace NaNs with the column medians
for i in range(adata.X.shape[1]):
    adata.X[:, i] = np.where(np.isnan(adata.X[:, i]), medians[i], adata.X[:, i])

adata.X = 1-adata.X
adata.X = sparse.csr_matrix(adata.X)
snap.tl.spectral(adata,features = None)
snap.pp.knn(adata)

def get_N_clusters(adata, n_cluster, range_min=0, range_max=3, max_steps=30, tolerance=0):
    """
    Tune the resolution parameter in clustering to make the number of clusters and the specified number as close as possible.
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    n_cluster
        Specified number of clusters.
    range_min
        Minimum clustering resolution for the binary search.
    range_max
        Maximum clustering resolution for the binary search.
    max_steps
        Maximum number of steps for the binary search.
    tolerance
        Tolerance of the difference between the number of clusters and the specified number.
    Returns
    -------
    adata
        AnnData object with clustering assignments in `adata.obs`:
    """ 
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        snap.tl.leiden(adata, resolution=this_resolution)
        print(this_resolution)
        print(adata.obs['leiden'])
        this_clusters = adata.obs['leiden'].nunique()
        if this_clusters > n_cluster+tolerance:
            this_max = this_resolution
        elif this_clusters < n_cluster-tolerance:
            this_min = this_resolution
        else:
            print(this_step)
            print("Succeed to find %d clusters at resolution %.3f."%(n_cluster, this_resolution))
            return adata
        this_step += 1
    print('Cannot find the number of clusters.')
    print(this_step)
    return adata

adata = get_N_clusters(adata, n_clusters)
sorted_cell_names = sorted(adata.obs_names)
cell2major = dict(zip(adata.obs.index, adata.obs["leiden"]))
major_type_list = [cell2major[cell] for cell in sorted_cell_names]
major_type_arr = np.array(major_type_list)
df_clustering = pd.DataFrame()
df_clustering['SnapATAC2-C'] = major_type_arr
df_clustering.to_csv(predcluster, header=True, index=False)

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
    print(f"ari=%.6f" % (ari))
    print(f"ami=%.6f" % (ami))
    print(f"nmi=%.6f" % (nmi))
    print(f"fmi=%.6f" % (fmi))
    print(f"comp=%.6f" % (comp))
    print(f"homo=%.6f" % (homo))
    df = pd.DataFrame({'index': ['ari', 'ami', 'nmi', 'fmi', 'comp', 'homo'], 'value': [ari, ami, nmi, fmi, comp, homo]})
    df.to_csv(savepath, sep='\t', header=True, index=False)

cluster_metrics(true_label, major_type_arr)
