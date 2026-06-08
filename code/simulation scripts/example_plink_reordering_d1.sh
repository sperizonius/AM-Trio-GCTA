## Example bash script of reordering GNAMES simulated data for design 1 into Trio-GCTA required data order with PLINK 1.9 and 2.0
## Sofie Perizonius

#!/bin/bash

iDESIGN=1

mapfile -t names < <(tr '\t' '\n' < "setting_list_DESIGN.${iDESIGN}.csv")

settinglist=()
for name in "${names[@]}"; do
  settinglist+=("$name")
done

for setting in "${settinglist[@]}"; do
  design_run_name="${setting}"
  children_name="children.${design_run_name}"
  parents_name="parents.${design_run_name}"
  merged_name="merged.${design_run_name}"
  reorder_name="reorder_${design_run_name}.txt"
  merged_reordered="merged_reordered.${design_run_name}"
  grm_name="grm.${design_run_name}"

  plink \
  --bfile "${children_name}" \
  --bmerge "${parents_name}" \
  --make-bed \
  --out "${merged_name}"

  plink2 \
  --bfile "${merged_name}" \
  --keep "${reorder_name}" \
  --indiv-sort file "${reorder_name}" \
  --make-bed \
  --out "${merged_reordered}"

  mv "${merged_reordered}"* to_analyze
done

mv "setting_list_DESIGN.${iDESIGN}.csv" to_analyze
mv "target_parameters_DESIGN.${iDESIGN}.csv" to_analyze
