## Script used to obtain subset of the original reorder file (used for reordering simulated datasets according to required Trio-GCTA 
## required order) that only contains IDs of the children for the subsequent GCTA analysis

import numpy as np
import pandas as pd

setting_list_d1 = pd.read_csv("simulated_data/DESIGN_1/to_analyze/setting_list_DESIGN.1.csv", sep='\t', header = None)
setting_list_d1 = setting_list_d1.T

for setting in setting_list_d1[0]:
    reorder_name =  'reorder_'+setting+'.txt'
    reorder_file = pd.read_csv('simulated_data/DESIGN_1/'+reorder_name, sep=r"\s+", header=None)
    reorder_subset = reorder_file.iloc[12000:18000, :]
    reorder_subset.to_csv('simulated_data/DESIGN_1/GCTA_subsets/'+'subset_children'+'reorder_.'+setting+'.txt', sep="\t", header=False, index=False)

setting_list_d2 = pd.read_csv("simulated_data/DESIGN_2/to_analyze/setting_list_DESIGN.2.csv", sep='\t', header = None)
setting_list_d2 = setting_list_d2.T

for setting in setting_list_d2[0]:
    reorder_name =  'reorder_'+setting+'.txt'
    reorder_file = pd.read_csv('simulated_data/DESIGN_2/'+reorder_name, sep=r"\s+", header=None)
    reorder_subset = reorder_file.iloc[12000:18000, :]
    reorder_subset.to_csv('simulated_data/DESIGN_2/GCTA_subsets/'+'subset_children'+'reorder_.'+setting+'.txt', sep="\t", header=False, index=False)


setting_list_d3 = pd.read_csv("simulated_data/DESIGN_3/to_analyze/setting_list_DESIGN.3.csv", sep='\t', header = None)
setting_list_d3 = setting_list_d3.T

for setting in setting_list_d3[0]:
    reorder_name =  'reorder_'+setting+'.txt'
    reorder_file = pd.read_csv('simulated_data/DESIGN_3/'+reorder_name, sep=r"\s+", header=None)
    reorder_subset = reorder_file.iloc[12000:18000, :]
    reorder_subset.to_csv('simulated_data/DESIGN_3/GCTA_subsets/'+'subset_children'+'reorder_.'+setting+'.txt', sep="\t", header=False, index=False)


setting_list_d4 = pd.read_csv("simulated_data/DESIGN_4/to_analyze/setting_list_DESIGN.4.csv", sep='\t', header = None)
setting_list_d4 = setting_list_d4.T

for setting in setting_list_d4[0]:
    reorder_name =  'reorder_'+setting+'.txt'
    reorder_file = pd.read_csv('simulated_data/DESIGN_4/'+reorder_name, sep=r"\s+", header=None)
    reorder_subset = reorder_file.iloc[12000:18000, :]
    reorder_subset.to_csv('simulated_data/DESIGN_4/GCTA_subsets/'+'subset_children'+'reorder_.'+setting+'.txt', sep="\t", header=False, index=False)


setting_list_d2_nocovar = pd.read_csv("simulated_data/DESIGN_2_NOCOVAR/to_analyze/setting_list_DESIGN.2.csv", sep='\t', header = None)
setting_list_d2_nocovar = setting_list_d2_nocovar.T

for setting in setting_list_d2_nocovar[0]:
    reorder_name =  'reorder_'+setting+'.txt'
    reorder_file = pd.read_csv('simulated_data/DESIGN_2_NOCOVAR/'+reorder_name, sep=r"\s+", header=None)
    reorder_subset = reorder_file.iloc[12000:18000, :]
    reorder_subset.to_csv('simulated_data/DESIGN_2_NOCOVAR/GCTA_subsets/'+'subset_children'+'reorder_.'+setting+'.txt', sep="\t", header=False, index=False)


setting_list_d3_nocovar = pd.read_csv("simulated_data/DESIGN_3_NOCOVAR/to_analyze/setting_list_DESIGN.3.csv", sep='\t', header = None)
setting_list_d3_nocovar = setting_list_d3_nocovar.T

for setting in setting_list_d3_nocovar[0]:
    reorder_name =  'reorder_'+setting+'.txt'
    reorder_file = pd.read_csv('simulated_data/DESIGN_3_NOCOVAR/'+reorder_name, sep=r"\s+", header=None)
    reorder_subset = reorder_file.iloc[12000:18000, :]
    reorder_subset.to_csv('simulated_data/DESIGN_3_NOCOVAR/GCTA_subsets/'+'subset_children'+'reorder_.'+setting+'.txt', sep="\t", header=False, index=False)
