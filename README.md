# AM-Trio-GCTA
This repository contains all code used to obtain the results of the simulation study and the empirical application described in the thesis _On the Importance of Assortative Mating Bias in Genetic Modelling: An Assessment and Extension of the Trio-GCTA Model_. This includes the adapted versions of the GNAMES tool (https://github.com/devlaming/gnames) (van Kippersluis et al., 2023) and the GREMLmodels package (https://github.com/espenmei/GREMLModels.jl) (Eilertsen et al., 2026) that were used for this thesis.

## Adapted GNAMES tool
For this project, data were simulated using the Genetic-Nurture and Assortative-Mating-Effects Simulator (GNAMES) tool (https://github.com/devlaming/gnames) (van Kippersluis et al., 2023). To meet the data requirements of the current project, some functions were modified or added. The adapted GNAMES tool enables separate specification of maternal and paternal indirect genetic effects, rather than a single combined parental effect. It also allows users to specify the covariance between direct and indirect genetic effects, and to output not only the PLINK files for the final generation but also those of their parents, which are required for Trio-GCTA analyses. Any functions that were adapted or added relative to the original GNAMES tool are explicitly indicated in the code (script is named /code/simulation study/GNAMES_adapted.py). In addition, explanatory comments have been added to many functions. 

## Adapted GREMLmodels package
Trio-GCTA models can be fitted with the GREMLmodels package (https://github.com/espenmei/GREMLModels.jl) (Eilertsen et al., 2026). We have adapted a few functions to enable fitting the new variance decomposition that comes with AM-Trio-GCTA. To use the adapted GREMLmodels version, download the code from the original package and replace the gremlmodels.jl sourcecode with the script that can be found in this repository (/code/GREMLmodels_adapted.jl).

## Simulation study
All (example) scripts used for the simulation study can be found in **/code**/**simulation study**.

### Data generation and analysis
The script example_simulation_d1.py contains example code to simulate a full design (in this case Design 1), with multiple runs and settings. The script also contains code to save the target results based on the last simulated generations. After simulating the data, the data has to be reordered in the order that the Trio-GCTA models require (mothers, fathers, children, sorted on family number). For this to be done, reorder files have to be made and this is also done in the example_simulation_d1.py script with the function prep_for_plink (function can be found in sorting_functions.py script). Afterwards, the actual reordering of the data can be done with the example_plink_reordering_d1.sh script. Be aware that both PLINK 1.9 and PLINK 2.0 (Chang et al., 2015) are required to run this script.

#### Trio-GCTA analysis
The script /trio-gcta/example_model_fitting_d1.py contains the code to fit the standard version of Trio-GCTA (mode = :standard) or the AM-Trio-GCTA version (mode = :am_version) to all simulated datasets within one design. The functions used in this script are in /trio-gcta/fitting_trio_functions.jl. The script /trio-gcta/pipeline_all_designs.sh contains the code to run the pipeline of simulating, reordering, and analyzing the data with the Trio-GCTA models for multiple designs.

#### Standard GCTA analysis
To also apply standard GCTA (Yang et al., 2011) to the simulated data, a subset of the PLINK files has to be made that only contains the data for the children. These scripts are in the folder /gcta. The scripts in this folder should be used in the following order: 1) make_reorder_susbsets.py, 2) plink_subsets.sh, 3) make_phen_file.py, 4) gcta_analysis.sh.

The GCTA software that was used was downloaded from: https://yanglab.westlake.edu.cn/software/gcta/#Overview. 

## Empirical application
All (example) scripts used for the empirical application can be found in **/code**/**empirical application**. To use the example scripts in this repository, start with a reorder file (txt) file that contains the IDs of individuals that remained after data cleaning (see thesis for data cleaning procedure) in the right Trio-GCTA-required order, and a .csv or .txt file that contains the same IDs (incl family ID) and the residualized and standardized phenotypes. The reorder file can be used to obtain reordered PLINK files for the selected individuals from general PLINK genotype files (script /trio-gcta/plink_reordering_and_selection.sh). We precomputed a GRM and saved this GRM after reordering the data to save computation time when (re)fitting the models. Afterwards, the script /trio-gcta/empirical_model_fitting.jl was used to fit classical Trio-GCTA or AM-Trio-GCTA to the reordered data and GRM. This script uses functions that are in /trio-gcta/empirical_trio_functions.jl. The difference with the functions in this file versus the functions used for the simulated data is that these functions take a precomputed GRM and take the phenotype information from a prespecified path instead of from the .fam file. 

After applying the Trio-GCTA models to the data, we also applied standard GCTA to the empirical data. To do this, we also had to select just the children from the data. A similar order as described in the simulation study was used and the scripts are in the /gcta folder. The scripts in this folder should be used in the following order: 1) empirical_make_reorder_susbsets.py, 2) empirical_plink_subsets.sh, 3) empirical_make_phen_file.py, 4) empirical_gcta_analysis.sh.

Be aware that all scripts contain references to folders, directories, or scripts that might be different when using the scripts in other settings.

## Author
- Sofie Perizonius, contact details: s.perizonius3@student.vu.nl or sofie.perizonius@gmail.com (after 1-09-2026)

## References
Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: Rising to the challenge of larger and richer datasets. GigaScience, 4 (1), Article 7. https://doi.org/10.1186/s13742-015-0047-8

Eilertsen, E. M. (2026). GREMLModels.jl - a julia package for fitting variance component models with relationship matrices (Version 0.1.4). https://github.com/espenmei/GREMLModels.jl

van Kippersluis, H., Biroli, P., Pereira, R. D., Galama, T. J., von Hinke, S., Meddens, S. F. W., Muslimova, D., Slob, E. A., de Vlaming, R., & Rietveld, C. A. (2023). Overcoming attenuation bias in regressions using polygenic indices. Nature Communications, 14 (1), Article 4473. https://doi.org/10.1038/s41467-023-40069-4

Yang, J., Lee, S. H., Goddard, M. E., & Visscher, P. M. (2011). GCTA: A tool for genome-wide complex trait analysis. American Journal of Human Genetics, 88 (1), 76–82. https://doi.org/10.1016/j.ajhg.2010.11.011

