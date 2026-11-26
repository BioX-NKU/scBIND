import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys

h5adpath = sys.argv[1]
transformpath = sys.argv[2]

adata = sc.read(h5adpath)
for i in range(adata.n_obs):
    df = pd.DataFrame()
    meth_idx = np.isnan(adata.X[i, :])
    meth_idx = np.where(meth_idx == False)[0]
    meth_status = adata.X[i, meth_idx]
    df['chr'] = [1] * len(meth_idx)
    df['pos'] = meth_idx
    df['Met'] = meth_status
    df['Total'] = [1] * len(meth_idx)
    df.to_csv(os.path.join(transformpath, f"{adata.obs_names[i]}.gz"), header=None, index=None, sep='\t', compression='gzip')

