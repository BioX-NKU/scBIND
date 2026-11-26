import pandas as pd
import scanpy as sc
import anndata as ad
import os
import sys

output_VMR_matrix = sys.argv[1]
groundtruth = sys.argv[2]

X = pd.read_csv(os.path.join(output_VMR_matrix, "pca_results_mtx.csv"), header=0, index_col=None, sep=',')
adata = ad.AnnData(X.iloc[:, :-1], var=pd.DataFrame(index=X.columns[:-1]), obs=pd.DataFrame(index=X.iloc[:, -1]))
cell_type = pd.read_csv(groundtruth, index_col=None, header=0, sep='\t')
adata.obs['cell_type'] = cell_type.iloc[:, 1].values
adata.write_h5ad("constructed.h5ad")
