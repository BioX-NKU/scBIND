SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1
set -e

h5admet="/path/to/your_met.h5ad"
h5adtot="/path/to/your_tot.h5ad"
threshold_cover="61"
threshold_min_mean="0.2"
threshold_max_mean="0.8"
groundtruth="/path/to/true_clone_membership.txt.gz"
n_clusters="6"
predcluster="./result/scMET.csv"
savepath="./result/Evaluation_scMET.csv"

conda run -n scMET python adjust.py $h5admet $h5adtot $threshold_cover $threshold_min_mean $threshold_max_mean > log 2>&1
conda run -n scMET Rscript run.R >> log 2>&1
conda run -n scMET python run.py $h5admet $h5adtot $groundtruth $n_clusters $predcluster $savepath >> log 2>&1

