SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1
set -e

input_path="/path/to/your/data-F/"
groundtruth="/path/to/your/true_clone_membership.txt.gz"
predcluster="./result/scMelody.csv"
savepath="./result/Evaluation_scMelody.csv"

conda run -n scMelody Rscript data_generation.r $input_path > log 2>&1
conda run -n scMelody python cluster.py $input_path $groundtruth $predcluster $savepath >> log 2>&1

