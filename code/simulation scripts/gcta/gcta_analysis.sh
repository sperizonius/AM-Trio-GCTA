## Script to perform GCTA analysis on children subsets of simulated datasets
## Sofie Perizonius - 9-6-2026

#!/bin/bash

for iDESIGN in 1 2 3 4; do
        mapfile -t names < <(tr '\t' '\n' < "simulated_data/DESIGN_${iDESIGN}/to_analyze/setting_list_DESIGN.${iDESIGN}.csv")

        settinglist=()
        for name in "${names[@]}"; do
        settinglist+=("$name")
        done

        plink_dir="simulated_data/DESIGN_${iDESIGN}/GCTA_subsets"
        grm_dir="simulated_data/DESIGN_${iDESIGN}/GCTA_subsets/grms"
        phenotype_dir="simulated_data/DESIGN_${iDESIGN}/GCTA_subsets/phenotypes"
        out_dir="simulation_results/DESIGN_${iDESIGN}/GCTA_results/h2_results"

        for setting in "${settinglist[@]}"; do
        design_run_name="${setting}"
        children_name="subset_children.${design_run_name}"
        phenotype_name="subset_children.${design_run_name}.phen"
        out_name="GRM_children.${design_run_name}"
        results_name="GCTA_results.${design_run_name}"

        ./gcta64 --bfile "${plink_dir}/${children_name}" \
                --make-grm \
                --pheno "${phenotype_dir}/${phenotype_name}" \
                --reml \
                --out "${grm_dir}/${out_name}" --thread-num 20

        ./gcta64 --grm  "${grm_dir}/${out_name}" \
                --pheno "${phenotype_dir}/${phenotype_name}" \
                --reml --out "${out_dir}/${results_name}" --thread-num 20
        done
done
