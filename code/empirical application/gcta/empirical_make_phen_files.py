## Script to make .phen files containing the phenotypes of the children for the GCTA analysis
## Sofie Perizonius - 9-6-2026

import numpy as np
import pandas as pd

def make_phen_children(fam_path, phenotype_path):    
    # read fam file 
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None)
    fam = pd.DataFrame(fam)

    phe_csv = pd.read_csv(phenotype_path)
    n_trait = int(phe_csv.shape[0] / 3)
    phe_column =  phe_csv.iloc[2*n_trait : 3*n_trait + 1, 3]
    phe = fam.loc[:, [0,1]]
    phe[2] = phe_column.values    

    return phe

phen_bw = make_phen_children('for_gcta/birth_weight_children.fam', "for_trio_gcta/birth_weight_trios.csv")
phen_bw.to_csv('for_gcta/birth_weight_children.phe', sep="\t", header=False, index=False)

phen_bl = make_phen_children('for_gcta/birth_length_children.fam', "for_trio_gcta/birth_length_trios.csv")
phen_bl.to_csv('for_gcta/birth_length_children.phe', sep="\t", header=False, index=False)

phen_height = make_phen_children('for_gcta/height_children.fam', "for_trio_gcta/height_trios.csv")
phen_height.to_csv('for_gcta/height_children.phe', sep="\t", header=False, index=False)

phen_edu = make_phen_children('for_gcta/edu_children.fam', "for_trio_gcta/edu_trios.csv")
phen_edu.to_csv('for_gcta/edu_children.phe', sep="\t", header=False, index=False)
