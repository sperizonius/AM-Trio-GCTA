## Script to make a subset of the original simulated datasets that were used in the Trio-GCTA analyses that only contains the children
## for the GCTA analysis

#!/bin/bash

# iDESIGN=1
# iDESIGN=2
# iDESIGN=3
# iDESIGN=4

# iDESIGN=2
# iDESIGN_nocovar=2_NOCOVAR

# iDESIGN=3
# iDESIGN_nocovar=3_NOCOVAR

mapfile -t names < <(tr '\t' '\n' < "simulated_data/DESIGN_${iDESIGN_nocovar}/to_analyze/setting_list_DESIGN.${iDESIGN}.csv")

settinglist=()
for name in "${names[@]}"; do
  settinglist+=("$name")
done

base_dir="simulated_data/DESIGN_${iDESIGN_nocovar}"
reorder_dir="simulated_data/DESIGN_${iDESIGN_nocovar}/GCTA_subsets"

for setting in "${settinglist[@]}"; do
    design_run_name="${setting}"
    children_name="children.${design_run_name}"
    reorder_name="subset_childrenreorder_.${design_run_name}.txt"
    out_name="subset_children.${design_run_name}"

    plink2 --bfile "${base_dir}/${children_name}" \
            --keep "${reorder_dir}/${reorder_name}" \
            --make-bed \
            --out "${reorder_dir}/${out_name}"
    done

