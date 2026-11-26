import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys

h5admet = sys.argv[1]
h5adtot = sys.argv[2]
threshold_cover = int(sys.argv[3])
threshold_min_mean = float(sys.argv[4])
threshold_max_mean = float(sys.argv[5])

adata_met = sc.read(h5admet)
adata_tot = sc.read(h5adtot)
gene_counts = (~np.isnan(adata_tot.X)).sum(axis=0)
gene_indices = np.where(gene_counts >= threshold_cover)[0] # scmet required features covering over 3 cells
adata_tot = adata_tot[:, gene_indices].copy()
adata_met = adata_met[:, gene_indices].copy()

# exclude features with extremely low variances, otherwise BB model cannot convergent
X = adata_met.X / adata_tot.X
means = np.nanmean(X, axis=0)
adata_tot = adata_tot[:, (means > threshold_min_mean) & (means < threshold_max_mean)].copy()
adata_met = adata_met[:, (means > threshold_min_mean) & (means < threshold_max_mean)].copy()

cols, rows = np.where(~np.isnan(adata_tot.X.T))
feature_names = adata_tot.var_names[cols]
cell_names = adata_tot.obs_names[rows]
vals_tot = np.asarray(adata_tot.X[rows, cols]).flatten()
vals_met = np.asarray(adata_met.X[rows, cols]).flatten()

df = pd.DataFrame({'Feature': feature_names, 'Cell': cell_names, 'total_reads': vals_tot, 'met_reads': vals_met})
df.to_csv('intermediate.csv', header=True, index=None, sep='\t')

