import sys
import scanpy as sc
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
import numpy as np

def calculate_cosine(x1, x2):
    mask = ~np.isnan(x1) & ~np.isnan(x2)
    x1_mask = x1[mask]
    x2_mask = x2[mask]
    if len(x1_mask) == 0:
        return np.nan
    return 1 - cosine(x1_mask, x2_mask)

def calculate_pearson(x1, x2):
    mask = ~np.isnan(x1) & ~np.isnan(x2)
    x1_mask = x1[mask]
    x2_mask = x2[mask]
    if len(x1_mask) == 0:
        return np.nan
    correlation, _ = pearsonr(x1_mask, x2_mask)
    return correlation

def calculate_hamming(x1, x2):
    mask = ~np.isnan(x1) & ~np.isnan(x2)
    x1_mask = x1[mask]
    x2_mask = x2[mask]
    if len(x1_mask) == 0:
        return np.nan
    return 1 - np.mean(x1_mask != x2_mask)

h5adpath = sys.argv[1]
adata = sc.read(h5adpath)
data = (adata.X >= 0.5).astype(int)

cosine_distance = np.zeros((data.shape[0], data.shape[0]))
pearson_distance = np.zeros((data.shape[0], data.shape[0]))
hamming_distance = np.zeros((data.shape[0], data.shape[0]))

for i in range(data.shape[0]):
    print(i)
    for j in range(data.shape[0]):
        cosine_distance[i][j] = calculate_cosine(data[i], data[j])
        pearson_distance[i][j] = calculate_pearson(data[i], data[j])
        hamming_distance[i][j] = calculate_hamming(data[i], data[j])

np.savetxt("./distance/dism_cosine.csv", cosine_distance, delimiter=',')
np.savetxt("./distance/dism_pearson.csv", pearson_distance, delimiter=',')
np.savetxt("./distance/dism_hamming.csv", hamming_distance, delimiter=',')

