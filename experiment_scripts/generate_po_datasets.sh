module load R
data_reps=96
#alpha="-2 -1.5 -1 -0.5"
alpha="-2"
beta="0.5"
gamma="1"
#gamma="-2 -1 -0.5 1"
#delta="0.25"
delta="-0.25 0.25 2"
# Might be nice to write job output to the experimental run dir but it's nice to leave that dir created by the R script
# so as to seperate the HPC and local run capabilities.
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
for g in $gamma
do
for d in $delta
do
Rscript $WORKDIR/generate_po_data.R --working_dir=$WORKDIR --exp_name=po_gen --alpha=$alpha --beta=$beta --gamma=$g --delta=$d --intensity_func="donut" --bias_fun="exponential" --data_reps=$data_reps
done
done
