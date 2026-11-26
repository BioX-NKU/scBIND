SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" || exit 1

h5adpath="/path/to/your.h5ad"
transformpath="/path/to/your/transformed_h5ad/"

conda run -n Melissa python build_Region.py $h5adpath $transformpath > log 2>&1
