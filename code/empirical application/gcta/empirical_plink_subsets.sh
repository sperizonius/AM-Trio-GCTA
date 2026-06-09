## Script to make subsets of the empirical data for only children with PLINK
## Sofie Perizonius - 9-9-2026

#!/bin/bash

plink2 --bfile "for_trio_gcta/birth_weight_trios" \
    --keep "for_trio_gcta/birth_weight_children.txt" \
    --make-bed \
    --out "for_gcta/birth_weight_children"

plink2 --bfile "for_trio_gcta/birth_length_trios" \
    --keep "for_trio_gcta/birth_length_children.txt" \
    --make-bed \
    --out "for_gcta/birth_length_children"

plink2 --bfile "for_trio_gcta/height_trios" \
    --keep "for_trio_gcta/height_children.txt" \
    --make-bed \
    --out "for_gcta/height_children"

plink2 --bfile "for_trio_gcta/edu_trios" \
    --keep "for_trio_gcta/edu_children.txt" \
    --make-bed \
    --out "for_gcta/edu_children"

