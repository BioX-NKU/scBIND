bed_path="/path/to/your/data/"
h5adpath="/path/to/your/save/"
chrom="0"
pos="1"
met="5"
tot="4"
species="HUMAN"
mm_refgenome="None"
hg_refgenome="hg19"
window_width="10000"
groundtruth="/path/to/your/true_clone_membership.txt.gz"

conda run -n scMET python build_scMET.py $bed_path $h5adpath $chrom $pos $met $tot $species $mm_refgenome $hg_refgenome $window_width $groundtruth > log 2>&1

