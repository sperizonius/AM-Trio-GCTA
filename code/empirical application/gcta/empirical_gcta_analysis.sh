## Script to compute GRMs for the subsets of the data containing only the children and to perform standard GCTA analyses on this data
## Sofie Perizonius - 9-6-2026

#!/bin/bash

# birth weight
./gcta64 --bfile "for_gcta/birth_weight_children" \
                --make-grm \
                --reml \
                --out "for_gcta/GRM_birth_weight_children" --thread-num 20

./gcta64 --grm  "for_gcta/GRM_birth_weight_children" \
                --pheno "for_gcta/birth_weight_children.phe" \
                --reml --out "for_gcta/results/birth_weight_GCTA_results" --thread-num 20

# birth length
./gcta64 --bfile "for_gcta/birth_length_children" \
                --make-grm \
                --reml \
                --out "for_gcta/GRM_birth_length_children" --thread-num 20

./gcta64 --grm  "for_gcta/GRM_birth_length_children" \
                --pheno "for_gcta/birth_length_children.phe" \
                --reml --out "for_gcta/results/birth_length_GCTA_results" --thread-num 20

# height
./gcta64 --bfile "for_gcta/height_children" \
                --make-grm \
                --reml \
                --out "for_gcta/GRM_height_children" --thread-num 20

./gcta64 --grm  "for_gcta/GRM_height_children" \
                --pheno "for_gcta/height_children.phe" \
                --reml --out "for_gcta/results/height_GCTA_results" --thread-num 20

# edu
./gcta64 --bfile "for_gcta/edu_children" \
                --make-grm \
                --reml \
                --out "for_gcta/GRM_edu_children" --thread-num 20

./gcta64 --grm  "for_gcta/GRM_edu_children" \
                --pheno "for_gcta/edu_children.phe" \
                --reml --out "for_gcta/results/edu_GCTA_results" --thread-num 20
