# Evaluation

All scripts used to evaluate the results of the optimal design algorithm are evailable in the ```evaluation``` directory within this repository. Typically, the worklflow of this code is designed to be run until the full design algorithm is completed. Then, the evaluation scripts can be run to generate tables and figures related to the set of experiments run for a particular set of parameters and experimental options. Furthermore, the optimal designs themselves can be found in the results directory within the saved files for a particular experiment. 

The instructions here describe how to recreate the figures and tables presented in the paper associated with this work. Beyond this, for your particular use-case and experimental settings, the scripts will need to be manually altered. 

## Scripts

For each experiment reported in the paper, there is a script to run the corresponding evaluation and generate the figures and tables shown in the paper. There are also a few other experimental evaluation scripts included in the code here that were not able to be reported in the paper due to lack of space.

```evaluation/read_job_output_exp_1.R``` will generate the results for the first experiment reported in the paper, in which including presence data was compared to results when it was not included and when the number of visits to each site was held constant.

EDIT from here: Say that you need the data, and such
```evaluation/read_job_output_exp_2.R```