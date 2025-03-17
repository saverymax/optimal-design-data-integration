#module load R
data_reps=96
alpha="-2"
beta="0.5 1"
gamma="1"
delta="0.25 2"
deviation=5
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
for d in $delta
do
for b in $beta
do
Rscript $WORKDIR/generate_po_data_misspec.R --working_dir=$WORKDIR --save_dir="misspec" --exp_name=po_gen --alpha=$alpha --beta=$b --gamma=$gamma --delta=$d --epsilon=1.5 --intensity_func="donut" --sd=$deviation --bias_fun="exponential" --data_reps=$data_reps --gamma_reps=200
done
done
