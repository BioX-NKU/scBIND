import pandas as pd
import numpy as np
import sys

groundtruth = sys.argv[1]
predcluster = sys.argv[2]
savepath = sys.argv[3]

true_label = pd.read_csv(groundtruth, sep='\t', header=0, index_col=None, compression='gzip').iloc[:, 1].values
print(true_label)
num_cells = len(true_label)

df = pd.read_csv(predcluster, sep='\t', header=0, index_col=None)
pred = df.iloc[:, 0].values
print(pred)

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

cluster_metrics(true_label, pred)
