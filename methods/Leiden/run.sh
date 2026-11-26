SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

h5adpath="/path/to/your.h5ad"
predcluster="./result/Leiden.csv"
savepath="./result/Evaluation_Leiden-C.csv"

conda run -n Louvain python Leiden.py $h5adpath $predcluster $savepath > log 2>&1

