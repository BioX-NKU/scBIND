SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1
set -e

h5adpath="/path/to/your.h5ad"
groundtruth="/path/to/true_clone_membership.txt.gz"
predcluster="./result/scMelodyRegion.csv"
savepath="./result/Evaluation_scMelodyRegion.csv"

conda run -n scMelody python data_generation.py $h5adpath > log 2>&1
conda run -n scMelody python cluster.py $h5adpath $groundtruth $predcluster $savepath >> log 2>&1

