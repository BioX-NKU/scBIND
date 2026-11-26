SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

h5adpath="/path/to/your.h5ad"
groundtruth="/path/to/your/true_clone_membership.txt.gz"
predcluster="./result/SnapATAC2-C.csv"
savepath="./result/Evaluation_SnapATAC2-C.csv"
n_clusters="6"

conda run -n SnapATAC2 python run.py $groundtruth $h5adpath $predcluster $savepath $n_clusters > log 2>&1

