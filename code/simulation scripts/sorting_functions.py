## Functions to sort .fam, .phe, .bim, and .bed files from GNAMES simulator into format that works with Trio-GCTA models
## Sofie Perizonius - 9-6-2026

import numpy as np
import pandas as pd

def phenotypes_to_fam(designname, children_or_parents = 'C'):
    """
    Takes last column from .phe file and pastes this into .fam file instead of original phenotype column in .fam file
    """
    # read phenotype file 
    if children_or_parents == 'C':
        phename = 'children.'+designname+'.phe'
        famname = 'children.'+designname+'.fam'
    
    if children_or_parents == 'P':
        phename = 'parents.'+designname+'.phe'
        famname = 'parents.'+designname+'.fam'
    
    phe = pd.read_csv(phename, sep=r"\s+", header=None)
    phe = pd.DataFrame(phe)

    # extract last column
    phenotypes = pd.Series(phe.iloc[:, -1])

    # read fam file 
    fam = pd.read_csv(famname, sep=r"\s+", header=None)
    fam = pd.DataFrame(fam)

    # safety check
    if len(phenotypes) != len(fam):
        raise ValueError("Phenotype file and Family file have different number of rows.")

    # Replace last column in fam file
    fam.loc[:, -1] = phenotypes
    fam = fam.drop([5], axis = 1)

    return fam

def add_FID(designname, n_children_per_fam = 1):
    """
    Creates family IDs for children and assigns same family IDs to parents of these children.
    """
    Cfamname = 'children.'+designname+'.fam'
    Pfamname = 'parents.'+designname+'.fam'

    children_fam = pd.read_csv(Cfamname, sep=r"\s+", header=None)
    parents_fam = pd.read_csv(Pfamname, sep=r"\s+", header=None)

    # children get a family number

    # first sort on parental ids such that children with the same parents appear together
    children_fam_sorted = children_fam.sort_values(by = [2, 3])

    n_fams = int(len(children_fam_sorted) / n_children_per_fam)
    family_ids = np.repeat(range(1, n_fams + 1), n_children_per_fam)
    children_fam_sorted.loc[:, 0] = family_ids

    # add family ids to children_fam such that original order of children_fam remains
    # merge on individual ID (column 1)
    children_fam = children_fam.merge(children_fam_sorted[[0, 1]],  on=1, how='left', suffixes=('', '_new'))

    # replace column 0 (zeros) with the merged family IDs
    children_fam.iloc[:,0] = children_fam['0_new']

    # crop helper column
    children_fam = children_fam.drop(columns=['0_new'])

    # create long format parent → family mapping
    father_df = children_fam.iloc[:, [0, 2]].rename(columns={2: 1})
    mother_df = children_fam.iloc[:, [0, 3]].rename(columns={3: 1})

    parent_map = pd.concat([father_df, mother_df])

    # remove duplicates (keep first occurrence per parent)
    parent_map = parent_map.drop_duplicates(subset=1)

    # merge with parents_fam on parental ID (column 1)
    parents_fam = parents_fam.merge(parent_map, on=1, how='left')

    # replace family ID (column 0) where match exists
    parents_fam.iloc[:,0] = parents_fam.iloc[:, 6]

    # drop helper column
    parents_fam = parents_fam.iloc[:, :-1]    

    # set parental ids in parents column to 0
    parents_fam.iloc[:, 2] = 0
    parents_fam.iloc[:, 3] = 0

    return children_fam, parents_fam

def make_reorder_file(designname, n_children_per_fam = 1):
    """
    Creates a "reorder" file that can be used in PLINK to reorder plink files in the way that they can be used in trio-GCTA.
    """
    Cfamname = 'children.'+designname+'.fam'
    Pfamname = 'parents.'+designname+'.fam'

    children_fam = pd.read_csv(Cfamname, sep=r"\s+", header=None)
    parents_fam = pd.read_csv(Pfamname, sep=r"\s+", header=None)

    # make reorder file based on family ids
    reorder_children = children_fam.iloc[:, 0:2]

    # randomly select one child per family such that the reorder file only contains one child per family if n_children_per_fam > 1
    if n_children_per_fam > 1:
            # randomly shuffle rows
            reorder_children = reorder_children.sample(frac=1, random_state=None)

            # keep only the first child per family (column 0 = family ID) (because first randomly shuffled, child that is picked is random)
            reorder_children = reorder_children.drop_duplicates(subset=0, keep='first').reset_index(drop=True)

    # sort reorder children set based on fam ID
    reorder_children = reorder_children.sort_values(by = 0)

    # sort reorder parents set based on fam ID and such that mothers come before fathers
    reorder_parents = parents_fam.iloc[:, [0, 1, 4]]
    reorder_parents = reorder_parents.sort_values(by = [4, 0], ascending = [False, True])
    reorder_parents = reorder_parents.iloc[:, [0, 1]]

    # remove individuals with family id 0 as they don't have children
    reorder_parents = reorder_parents[reorder_parents.iloc[:, 0] != 0]

    # create final file that contains the order that is needed for trio-GCTA, with index 1-n corresponding to mother, n-2n corresponding to fathers, and 2n-3n corresponding to children
    reorder = pd.concat([reorder_parents, reorder_children])
    return reorder

def prep_for_plink(dName, n_children_per_fam = 1):
    """
    Sequence of functions to add phenotypes to fam files, add family IDs to children and parent files and make reorder file for PLINK.
    """
    children_famname = 'children.'+ dName +'.fam'
    parents_famname = 'parents.'+ dName +'.fam'

    # Step 1: copy original fam files of parents and children and rename them with prefix original_
    original_children_fam = pd.read_csv(children_famname, sep=r"\s+", header=None)
    original_parents_fam = pd.read_csv(parents_famname, sep=r"\s+", header=None)
    original_children_fam.to_csv('original_' + children_famname, sep="\t", header=False, index=False)
    original_parents_fam.to_csv('original_'+ parents_famname, sep="\t", header=False, index=False)

    # Step 2: add phenotypes from .phe files to .fam files
    Cfam_file_incl_phe = phenotypes_to_fam(dName, children_or_parents = 'C')
    Cfam_file_incl_phe.to_csv(children_famname, sep="\t", header=False, index=False)

    Pfam_file_incl_phe = phenotypes_to_fam(dName, children_or_parents = 'P')
    Pfam_file_incl_phe.to_csv(parents_famname, sep="\t", header=False, index=False)

    # Step 3: add family IDs to PLINK files 
    Cfam_file_incl_phe_ID, Pfam_file_incl_phe_ID = add_FID(dName, n_children_per_fam)
    Cfam_file_incl_phe_ID.to_csv('children.'+ dName +'.fam', sep="\t", header=False, index=False)
    Pfam_file_incl_phe_ID.to_csv('parents.'+ dName +'.fam', sep="\t", header=False, index=False)

    # Step 4: create reorder file that can be used to obtain trio-GCTA required order
    reorder_file = make_reorder_file(dName, n_children_per_fam)
    reorder_file.to_csv('reorder_'+ dName +'.txt', sep="\t", header=False, index=False)


