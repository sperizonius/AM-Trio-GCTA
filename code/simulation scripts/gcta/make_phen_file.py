## Script to make .phen files containing just the phenotypes of the children for the GCTA analysis

import numpy as np
import pandas as pd

def make_phen(designname, path):
    # read phenotype file 
    famname = 'subset_children.'+designname+'.fam'
    
    # read fam file 
    fam = pd.read_csv(path + famname, sep=r"\s+", header=None)
    fam = pd.DataFrame(fam)

    phe = fam.iloc[:, [0,1,5]]

    return phe

setting_list_d1 = pd.read_csv("simulated_data/DESIGN_1/to_analyze/setting_list_DESIGN.1.csv", sep='\t', header = None)
setting_list_d1 = setting_list_d1.T

for setting in setting_list_d1[0]:
    new_phen_file = make_phen(setting, 'simulated_data/DESIGN_1/GCTA_subsets/')
    new_phen_file.to_csv('simulated_data/DESIGN_1/GCTA_subsets/phenotypes/'+'subset_children.'+setting+'.phen', sep="\t", header=False, index=False)

setting_list_d2 = pd.read_csv("simulated_data/DESIGN_2/to_analyze/setting_list_DESIGN.2.csv", sep='\t', header = None)
setting_list_d2 = setting_list_d2.T

for setting in setting_list_d2[0]:
    new_phen_file = make_phen(setting, 'simulated_data/DESIGN_2/GCTA_subsets/')
    new_phen_file.to_csv('simulated_data/DESIGN_2/GCTA_subsets/phenotypes/'+'subset_children.'+setting+'.phen', sep="\t", header=False, index=False)

setting_list_d3 = pd.read_csv("simulated_data/DESIGN_3/to_analyze/setting_list_DESIGN.3.csv", sep='\t', header = None)
setting_list_d3 = setting_list_d3.T

for setting in setting_list_d3[0]:
    new_phen_file = make_phen(setting, 'simulated_data/DESIGN_3/GCTA_subsets/')
    new_phen_file.to_csv('simulated_data/DESIGN_3/GCTA_subsets/phenotypes/'+'subset_children.'+setting+'.phen', sep="\t", header=False, index=False)

setting_list_d4 = pd.read_csv("simulated_data/DESIGN_4/to_analyze/setting_list_DESIGN.4.csv", sep='\t', header = None)
setting_list_d4 = setting_list_d4.T

for setting in setting_list_d4[0]:
    new_phen_file = make_phen(setting, 'simulated_data/DESIGN_4/GCTA_subsets/')
    new_phen_file.to_csv('simulated_data/DESIGN_4/GCTA_subsets/phenotypes/'+'subset_children.'+setting+'.phen', sep="\t", header=False, index=False)

setting_list_d2_nocovar = pd.read_csv("simulated_data/DESIGN_2_NOCOVAR/to_analyze/setting_list_DESIGN.2.csv", sep='\t', header = None)
setting_list_d2_nocovar = setting_list_d2_nocovar.T

for setting in setting_list_d2_nocovar[0]:
    new_phen_file = make_phen(setting, 'simulated_data/DESIGN_2_NOCOVAR/GCTA_subsets/')
    new_phen_file.to_csv('simulated_data/DESIGN_2_NOCOVAR/GCTA_subsets/phenotypes/'+'subset_children.'+setting+'.phen', sep="\t", header=False, index=False)

setting_list_d3_nocovar = pd.read_csv("simulated_data/DESIGN_3_NOCOVAR/to_analyze/setting_list_DESIGN.3.csv", sep='\t', header = None)
setting_list_d3_nocovar = setting_list_d3_nocovar.T

for setting in setting_list_d3_nocovar[0]:
    new_phen_file = make_phen(setting, 'simulated_data/DESIGN_3_NOCOVAR/GCTA_subsets/')
    new_phen_file.to_csv('simulated_data/DESIGN_3_NOCOVAR/GCTA_subsets/phenotypes/'+'subset_children.'+setting+'.phen', sep="\t", header=False, index=False)
