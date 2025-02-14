module load R
data_reps=96
alpha="-2"
beta="0.5 1"
gamma="1"
delta="2"
deviation=5
# Might be nice to write job output to the experimental run dir but it's nice to leave that dir created by the R script
# so as to seperate the HPC and local run capabilities.
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
for d in $delta
do
for b in $beta
do
Rscript $WORKDIR/generate_po_data_misspec.R --working_dir=$WORKDIR --exp_name=po_gen_peak --alpha=$alpha --beta=$b --gamma=$g --delta=$d --intensity_func="donut" --sd=$deviation --bias_fun="exponential" --data_reps=$data_reps
done
done
done
