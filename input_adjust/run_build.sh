bed_path="/path/to/your/data/"
h5adpath="/path/to/your.h5ad"
chrom="0"
pos="1"
met="4"
tot="6"
species="MOUSE"
mm_refgenome="mm10"
hg_refgenome="None"
window_width="10000"
groundtruth="/path/to/your/true_clone_membership.txt.gz"

conda run -n Louvain python build.py $bed_path $h5adpath $chrom $pos $met $tot $species $mm_refgenome $hg_refgenome $window_width $groundtruth > log 2>&1

