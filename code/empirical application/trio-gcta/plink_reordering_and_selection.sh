## Script to filter and reorder the data from the individuals that should be analyzed within the Trio-GCTA models
## Sofie Perizonius - 9-6-2026

#!/bin/bash
base_dir="/data/sofie/empirical_trio/GENOTYPE_RUNNING_FILES"
out_dir="/data/sofie/empirical_trio/for_trio_gcta"

cd /data/sofie/empirical_trio

# filter and sort trios for final trio gcta analysis
plink2 --bfile "${base_dir}/merged_cleaned_pot_trios" \
    --keep "${out_dir}/reorder_birth_length.txt" \
    --indiv-sort file "${out_dir}/reorder_birth_length.txt" \
    --make-bed \
    --out "${out_dir}/birth_length_trios"

plink2 --bfile "${base_dir}/merged_cleaned_pot_trios" \
    --keep "${out_dir}/reorder_birth_weight.txt" \
    --indiv-sort file "${out_dir}/reorder_birth_weight.txt" \
    --make-bed \
    --out "${out_dir}/birth_weight_trios"

plink2 --bfile "${base_dir}/merged_cleaned_pot_trios" \
    --keep "${out_dir}/reorder_height.txt" \
    --indiv-sort file "${out_dir}/reorder_height.txt" \
    --make-bed \
    --out "${out_dir}/height_trios"

plink2 --bfile "${base_dir}/merged_cleaned_pot_trios" \
    --keep "${out_dir}/reorder_edu.txt" \
    --indiv-sort file "${out_dir}/reorder_edu.txt" \
    --make-bed \
    --out "${out_dir}/edu_trios"

