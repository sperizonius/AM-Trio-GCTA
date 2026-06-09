# AM-Trio-GCTA
This repository contains all code used to obtain the results of the simulation study and the empirical application described in the thesis _On the importance of assortative mating bias in genetic modelling: An assessment and extension of the Trio-GCTA model_. This includes the adapted versions of the GNAMES tool (https://github.com/devlaming/gnames) and the GREMLmodels package (https://github.com/espenmei/GREMLModels.jl) that were used for this thesis.

## Simulation study
### GNAMES tool
For this project, data were simulated using the Genetic-Nurture and Assortative-Mating-Effects Simulator (GNAMES) tool (van Kippersluis et al., 2023). To meet the data requirements of the current project, some functions were modified or added. The adapted GNAMES tool enables separate specification of maternal and paternal indirect genetic effects, rather than a single combined parental effect. It also allows users to specify the covariance between direct and indirect genetic effects, and to output not only the PLINK files for the final generation but also those of their parents, which are required for Trio-GCTA analyses. Any functions that were adapted or added relative to the original GNAMES tool are explicitly indicated in the code. In addition, explanatory comments have been added to many functions.

> beschrijf of er nog iets nodig is om te installeren etc. > dependencies

> beschrijf naam script!

### Cleaning and preprocessing

> beschrijf volgorde scripts

### GREMLmodels package
> beschrijf of er nog iets nodig is om te installeren etc. > dependencies

### Standard GCTA application on GNAMES data
> beschrijf of er nog iets nodig is om te installeren etc. > dependencies

> script 1: make_reorder_subsets.py > ## Script used to obtain subset of the original reorder file (used for reordering simulated datasets according to required Trio-GCTA required order) that only contains IDs of the children for the subsequent GCTA analysis
> script 2: plink_subsets.sh > ## Script to make a subset of the original simulated datasets that were used in the Trio-GCTA analyses that only contains the children for the GCTA analysis
> script 3: make_phen_file.py > ## Script to make .phen files containing just the phenotypes of the children for the GCTA analysis
> script 4: gcta_analysis.sh > ## Script to perform GCTA analysis on children subsets of simulated datasets


## Empirical application
### Cleaning and preprocessing

> script 1: data_cleaning.R > contains all steps taken to clean the data prior to analysis with GCTA/Trio-GCTA. Output contains a reorder file that can be used by plink with all individuals in the right order for Trio-GCTA that remained after datacleaning and csv files containing the phenotype information of these individuals in the same order

> script ?: empirical_trio_functions.jl > ## Script that contains functions to fit trio-gcta models, either in the "classical" way or the AM-adjusted version, on simulated data from GNAMES and functions to fit all simulated datasets of one design in one go

Difference with functions used on simulated data is that these functions that a precomputed grm and take the (residualized) phenotype information 
from a csv file

> before using any of the scripts, please check if the references to folder, working directories or other scripts are correct.

### Structure
This list contains the structure of this repository and what different folders or files contain:
- **/code**: 
    - **simulation study**:
    - **empirical application**:

## Author
- Sofie Perizonius, contact details: s.perizonius3@student.vu.nl or sofie.perizonius@gmail.com (after 1-09-2026)


## References
van Kippersluis, H., Biroli, P., Pereira, R. D., Galama, T. J., von Hinke, S., Meddens, S. F. W., Muslimova, D., Slob, E. A., de Vlaming, R., & Rietveld, C. A. (2023). Overcoming attenuation bias in regressions using polygenic indices. Nature Communications, 14 (1), Article 4473. https://doi.org/10.1038/s41467-023-40069-4
