import pandas as pd
import sys
import numpy as np
import pandas as pd

predcluster = sys.argv[1]

df = pd.read_csv(predcluster, sep='\t', header=0, index_col=0)
num_cells = len(df)
df = df.iloc[0:num_cells, :].values
pred = [None] * num_cells
for i in range(0, num_cells):
    Max = 0
    Maxj = 0
    for j in range(0, df.shape[1]):
        if df[i][j] > Max:
            Max = df[i][j]
            Maxj = j
    pred[i] = Maxj
pred = np.array(pred)
df_clustering = pd.DataFrame()
df_clustering['Melissa'] = pred
df_clustering.to_csv(predcluster, header=True, index=False)

