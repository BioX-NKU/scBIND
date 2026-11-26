SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1
set -e

input_path="/path/to/your/data-S/*.gz"
output_path="/path/to/your/compact_data"
output_VMRs="/path/to/your/VMR/VMRs.bed"
output_VMR_matrix="/path/to/your/VMR/VMR_matrix"
groundtruth="/path/to/your/true_clone_membership.txt.gz"
predcluster="./result/MethSCAn-C.csv"
savepath="./result/Evaluation_MethSCAn-C.csv"
n_clusters="6"

conda run -n MethSCAn methscan prepare $input_path $output_path > log 2>&1
conda run -n MethSCAn methscan smooth $output_path >> log 2>&1
conda run -n MethSCAn methscan scan --threads 16 $output_path $output_VMRs >> log 2>&1
conda run -n MethSCAn methscan matrix --threads 16 $output_VMRs $output_path $output_VMR_matrix >> log 2>&1
conda run -n MethSCAn Rscript get_result.r $output_VMR_matrix >> log 2>&1 
conda run -n MethSCAn python build.py $output_VMR_matrix $groundtruth >> log 2>&1
conda run -n MethSCAn python Leiden-C.py $predcluster $savepath $n_clusters >> log 2>&1

