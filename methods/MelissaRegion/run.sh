SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1
set -e

transformpath="/path/to/your/transformed_h5ad/"
groundtruth="/path/to/your/true_clone_membership.txt.gz"
predcluster="./result/MelissaRegion.csv"
savepath="./result/Evaluation_MelissaRegion.csv"
n_clusters="6"

conda run -n Melissa Rscript run_melissa.R $transformpath $predcluster $n_clusters > log 2>&1
conda run -n Melissa python adjust_result.py $predcluster >> log 2>&1
conda run -n Melissa python Performance_Evaluate.py $groundtruth $predcluster $savepath >> log 2>&1
