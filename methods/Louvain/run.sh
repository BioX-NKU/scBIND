SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

h5adpath="/path/to/your.h5ad"
predcluster="./result/Louvain.csv"
savepath="./result/Evaluation_Louvain.csv"

conda run -n Louvain python Louvain.py $h5adpath $predcluster $savepath > log 2>&1

