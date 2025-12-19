# Introduction
Welcome to the documentation for the code repository "optimal-design-data-integration". 

The aim of this project is to explore the integration of presence-absence and presence-only data when designing surveys for biodiversity monitoring purposes. In this documentation you can find details regarding the usage of the code for designing surveys, implementation of the optimal design algorithms, and data relevant for the project.

To get started, please clone or copy the repository and follow along with the documentation here. If you wish to process the data yourself, make sure the data provided in the figshare repository <https://doi.org/10.6084/m9.figshare.25450108> is unzipped. All source code is available in the ```optimal-design-data-integration``` directory. After following the installation instructions, please go to the data section and follow the setup there. 

**Important:** This project is written in R. However, it is written as a Command Line application. Running the code in Rstudio will not work unless the scripts are substantially modified. Users must use the ```Rscript``` command from the CLI to run the modules, with the correct arguments as documented in the following sections. On Windows it may be necessary to provide the path to the ```Rscript``` executable, such as, ```"C:/Programs/R/R-42~1.2/bin/x64/Rscript.exe"```.
