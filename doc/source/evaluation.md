# Evaluation

All scripts used to evaluate the results of the optimal design algorithm are evailable in the ```evaluation``` directory within this repository. The workflow of the evaluation code is set up to be run after the design algorithm is completed, across different settings of the parameters or algorithm. Then, the evaluation scripts can be run to generate tables and figures for the set of experiments across a particular set of parameters and options. 

The instructions here describe how to recreate the figures and tables presented in the paper associated with this work. Beyond this, for your particular use-case and experimental settings, the scripts will need to be manually altered if you would like to compare your own runs.

## Scripts

For each experiment reported in the paper, there is a script to run the corresponding evaluation and generate the figures and tables shown in the paper. 

For example, ```Rscript read_job_output_exp_1.R --base_dir=your/path/to/experimental_runs``` will generate the results for the first experiment reported in the paper, in which results were compared when presence-only data was and was not included and when the number of visits to each site was held constant. To run this script, you need to have run the design algorithm with all permutations presented in the paper. The script to do this is in the ```experiment_scripts``` directory, as explained in the [experiments section](experiments) of this documentation. **However**, we have provided all necessary result files to run the evaluation scripts here. They are structured in the ```experimental_runs``` directory under the name of each experiment. The subdirectories there contain ```.xlsx``` files for each run reported in the paper. The evalution scripts will generate the tables and figures from the paper using these files. However, it is also possible to repeat the experiments yourself using the scripts in ```experiment_scripts```, given you have enough compute. 

By navigating to the evalution directory in a Linux system and running
```Rscript read_job_output_exp_1.R --base_dir=your/path/to/experimental_runs > exp_1_table.txt```
a Latex table will be created and can be found at the end of the .txt file that matches Table 1 from the project's paper. Make sure that the base directory you provide really leads to the ```experimental_runs``` directory containing the ```exp_1``` directory. As stated above, within this ```exp_1``` directory are many subdirectories, each with the results from a particular experiment provided in ```results.xlsx```. It is these results from each experiment that the script uses to create the Latex table. 

The rest of the scripts can be run as follows:

```Rscript read_job_output_exp_2.R --base_dir=your/path/to/experimental_runs > exp_2_table.txt```
This will generate a table comparing models with and without PO data, when visits are and are not allowed to vary, over a few parameter setings of beta and delta.
A table comparing the total number of visits for each run is also generated.

```Rscript read_job_output_exp_3.R --base_dir=your/path/to/experimental_runs > exp_3_table.txt```
This will generate a table and a figure showing the effect of increasing the maximum number of visits to each site.

```Rscript read_job_output_exp_4.R --base_dir=your/path/to/experimental_runs > exp_4_table.txt```
This generates a table and figure showing the effect of increasing the number of total sites used in the optimal design.

```Rscript read_job_output_application_1.R --base_dir=your/path/to/experimental_runs > application_1_table.txt```
This generates a table comparing the optimal designs as the amount of PO data used in the site occupancy model is increased.

```Rscript read_job_output_application_2.R --base_dir=your/path/to/experimental_runs > application_2_table.txt```
Finally, this script generates a table comparing the effect of increasing visits and sites on the designs created using a site occupancy model with no PO data.

We have provided a bash script to run all of these in the evaluation directory, ```read_all_results.sh```. Give the full path to the experimental_runs directory provided with this code that contains all the experiments with the .xlsx files, and run the bash script. This will generate all tables and figures.

In the evaluation folder, we also provide the output of each script in .txt files, as presented in the paper. 

### Misspecification analysis
