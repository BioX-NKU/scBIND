import pandas as pd
import numpy as np
from scipy.stats import entropy

methods = ['Louvain-C', 'Leiden-C', 'EpiclomalBasic', 'EpiclomalRegion',
           'Melissa', 'MelissaRegion', 'scMelody-C', 'scMelodyRegion-C', 
           'scMET', 'MethSCAn-C', 'SnapATAC2-C']
n_methods = len(methods)
n_cell = 122
n_clusters = 6
missing_rate = 0.240576744

df = pd.DataFrame(np.zeros((n_cell, n_methods)), columns=methods)

for i in range(n_methods):
    if methods[i] == 'EpiclomalBasic':
        df_read = pd.read_csv(f'/path/to/your/result/Epiclomal.csv', header=0, index_col=None)
        pred = df_read.iloc[:, 0].values
    elif methods[i] == 'EpiclomalRegion':
        df_read = pd.read_csv(f'/path/to/your/result/Epiclomal.csv', header=0, index_col=None)
        pred = df_read.iloc[:, 1].values
    else:
        df_read = pd.read_csv(f'/path/to/your/result/{methods[i]}.csv', header=0, index_col=None)
        pred = df_read.iloc[:, 0].values
    df[methods[i]] = pred

print(df)

from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, fcluster

def negmm_consensus(df, tau=0.5, n_clusters=3, linkage_method='average'):
    n_samples = len(df)
    n_base = df.shape[1]
    
    # co-association matrix construction
    co_matrix = np.zeros((n_samples, n_samples))
    for col in df.columns:
        labels = df[col].values
        adj_matrix = (labels.reshape(-1, 1) == labels.reshape(1, -1)).astype(int)
        co_matrix += adj_matrix
    
    co_matrix = co_matrix / n_base
    co_matrix_filtered = co_matrix.copy()
    co_matrix_filtered[co_matrix_filtered < tau] = 0 # remove weak association
    distance_matrix = 1 - co_matrix_filtered
    
    # apply allomerative clustering
    consensus_clustering = AgglomerativeClustering(
        n_clusters=n_clusters,
        metric='precomputed',
        linkage=linkage_method
    )
    np.fill_diagonal(distance_matrix, 0)
    consensus_labels = consensus_clustering.fit_predict(distance_matrix)
    return consensus_labels, co_matrix_filtered

def iterative_negmm(df, initial_tau=0, tau_step=0.05, n_clusters=3):
    current_df = df.copy()
    tau = initial_tau
    
    while tau <= 0.96:
        print(f"\n=== tau={tau:.2f} ===")
        labels, co_mat = negmm_consensus(current_df, tau=tau, n_clusters=n_clusters)
        non_zero_ratio = np.count_nonzero(co_mat) / co_mat.size
        print(f"Non zero proportion: {non_zero_ratio:.2%}")
        
        final_labels = labels
        
        if non_zero_ratio < 0.4 * (1 - missing_rate): # stop weak association removal if associations are few
            print("Too few associations. Stop removals.")
            break
        
        tau += tau_step
    
    return final_labels

labels = iterative_negmm(df, n_clusters=n_clusters)
pd.DataFrame({'scBIND': labels}).to_csv(f'./result/scBIND.csv', header=True, index=False)
true_label = pd.read_csv(f'/path/to/your/true_clone_membership.txt.gz', header=0, index_col=None, compression='gzip', sep='\t').iloc[:, 1].values

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
    df.to_csv(f'./result/Evaluation_scBIND.csv', sep='\t', header=True, index=False)

cluster_metrics(true_label, labels)
