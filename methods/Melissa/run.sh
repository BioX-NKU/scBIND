SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1
set -e

input_path="/path/to/your/data-F/"
groundtruth="/path/to/your/true_clone_membership.txt.gz"
predcluster="./result/Melissa.csv"
savepath="./result/Evaluation_Melissa.csv"
n_clusters="6"

conda run -n Melissa Rscript run_melissa.R $input_path $predcluster $n_clusters > log 2>&1
conda run -n Melissa python adjust_result.py $predcluster >> log 2>&1
conda run -n Melissa python Performance_Evaluate.py $groundtruth $predcluster $savepath >> log 2>&1
