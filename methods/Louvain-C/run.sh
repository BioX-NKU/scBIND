SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

h5adpath="/path/to/your.h5ad"
predcluster="./result/Louvain-C.csv"
savepath="./result/Evaluation_Louvain-C.csv"
n_clusters="6"

conda run -n Louvain python Louvain-C.py $h5adpath $predcluster $savepath $n_clusters > log 2>&1

