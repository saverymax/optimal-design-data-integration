# Experiments

This section briefly describes how to run the experiments showed in the paper associated with this work. Each experiment contains runs of the of optimal design algorithm for different parameter settings and experimental options. To run these, it is necessary to have access to an HPC cluster or otherwise large machine, as the optimal design algorithm uses intensive parallel computing and you will need to run this over many different experimental settings. That being said, we have provided the .xlsx files containing the results of the experiments presented in the paper, stored in the ```experimental_runs``` directory. Thus the instructions in this section are not necessary if you are interested in just rerunning the evaluation and recreating the tables and figures shown in the paper. In that case, the [evaluation section](evaluation) describes the steps to do so.

Otherwise, navigate to the```experiment_scripts``` directory. There are 6 bash scripts to generate the jobs for each experiment. These create ```.pbs``` files for running on Torque-based scheduler system. If your HPC uses slurm, these will not work, but it is not difficult to modify them. The bash scripts are named
```
make_exp_1_param_perm.sh
make_exp_2_model_1-3.sh
make_exp_3_model_survey_effort.sh
make_exp_4_m_increase.sh
make_application_script_1.sh
make_application_script_2.sh
```
These will generate all the necessary job files in respective directories. The scripts will also create a bash file to run the .pbs jobs. These will be named
```
run_exp_1.sh
run_exp_2.sh
run_exp_3.sh
run_exp_4.sh
run_applied_model_comparison_1.sh
run_applied_model_comparison_2.sh
```
If you run any of these files, the optimal designs will be saved within the experimental_runs directory, or the directory that you specify within the make scripts. It is then possible to compare the results between different parameter and experimental options using the evaluation scripts described in the next section.


