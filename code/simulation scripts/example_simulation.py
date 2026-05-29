import numpy as np
import pandas as pd

from GNAMES_adapted import gnames
from sorting_functions import *

## DESIGN 1: No indirect genetic effects
iDESIGN=1

# settings
iN = 12000
iM = 5000
T = np.array((1, 3, 20), int) 
iC = 2
dHsqY  = 0.6
dPropGN_M = 0
dPropGN_F = 0
dCov_G_GNM = 0
dCov_G_GNF = 0
dCov_GNM_GNF = 0
RhoAM = np.array((0, 0.25, 0.5))

d = {"am": [0, 0.25, 0.5, 0, 0.25, 0.5, 0, 0.25, 0.5], "generations": [1, 1, 1, 3, 3, 3, 20, 20, 20]}
settings = pd.DataFrame(d)
iSETTINGS = settings.shape[0]

# number of runs
nRun=10

def make_seed(design, setting, run):
    # convert to zero-based indexing
    design  -= 1
    setting -= 1
    run     -= 1

    seed = (
        design |
        (setting << 2) |
        (run << 6) 
    )
    
    return seed

def decode_seed(seed):
    design  = (seed & 0b11) + 1          # 2 bits
    setting = ((seed >> 2) & 0b1111) + 1 # 4 bits
    run     = ((seed >> 6) & 0b1111111) + 1 # 7 bits

    return {
        "design": design,
        "setting": setting,
        "run": run
    }

# dataframe to collect outcome parameters from simulations
target_results = pd.DataFrame()
setting_list = np.array([])

for iSetting in range(iSETTINGS):
    for run in range(1, nRun + 1):
        print(f"Run: {run}")

        #get unique seed for given combination of run, setting, and design
        #bits 1–4  : design
        # bits 5–8  : setting
        # bits 9–16 : run
        RhoAM = settings.iloc[iSetting, 0]
        T = settings.iloc[iSetting, 1]

        iThisSeed=make_seed(iDESIGN, iSetting + 1, run)

        sName='DESIGN.'+str(iDESIGN)+'.RhoAM.'+'{:.2f}'.format(RhoAM)+'.iT.'+'{:}'.format(T)+'.RUN.'+str(run)
        setting_list = np.append(setting_list, sName)

        simulator=gnames(iN,iM,iC,dHsqY,dPropGN_M, dPropGN_F, \
                        dCov_G_GNM, dCov_G_GNF, dCov_GNM_GNF, \
                            dRhoAM=RhoAM,iSeed=iThisSeed, bRescale=False)
        simulator.Simulate(int(T))

        cName = 'children.'+ sName
        pName = 'parents.'+ sName
        simulator.MakeBed(cName)
        simulator.MakeBedParents(pName)

        # prepare data for plink by reordering etc.
        prep_for_plink(sName, n_children_per_fam = iC)

        # get results of simulation
        # Flatten all offspring vectors
        # Number of siblings per family
        iS = simulator.mGY.shape[0]

        # Broadcast maternal/paternal effects to each sibling
        vGN_M_full = np.tile(simulator.vGN_M, (iS, 1)) 
        vGN_F_full = np.tile(simulator.vGN_F, (iS, 1))  

        vGY    = simulator.mGY.flatten()      # direct genetic effects
        vGN_M_full_flat = vGN_M_full.flatten()
        vGN_F_full_flat = vGN_F_full.flatten()
        vY     = simulator.mY.flatten()       # phenotype
        vE = simulator.mEY.flatten()

        # Compute variances
        var_G    = np.var(vGY)    # direct genetic variance
        var_GNM  = np.var(vGN_M_full_flat)  # maternal GN variance
        var_GNF  = np.var(vGN_F_full_flat)  # paternal GN variance
        var_E    = np.var(vE)     # environmental variance

        # Compute covariances (population covariances)
        cov_G_GNM = np.cov(vGY, vGN_M_full_flat, bias=True)[0,1]
        cov_G_GNF = np.cov(vGY, vGN_F_full_flat, bias=True)[0,1]
        cov_GNM_GNF = np.cov(vGN_M_full_flat, vGN_F_full_flat, bias=True)[0,1]

        # Actual phenotype variance
        V_Y_actual = np.var(vY)

        # save target results to dataframe
        target_parameters_this_run = np.array([var_GNM, var_GNF, var_G, cov_GNM_GNF, cov_G_GNM, cov_G_GNF, var_E, V_Y_actual])
        target_results.loc[:, sName] = target_parameters_this_run

target_results = target_results.rename(index={0: "Maternal indirect effect",
        1: "Paternal indirect effect",
        2: "Offspring direct effect",
        3: "Maternal-paternal covariance",
        4: "Maternal-offspring covariance",
        5: "Paternal-offspring covariance",
        6: "Residual variance",
        7: "Phenotypic variance"})   

target_results.to_csv('target_parameters_'+'DESIGN.'+str(iDESIGN) + '.csv', sep="\t", header=True, index=True)
setting_list.tofile('setting_list_' + 'DESIGN.' + str(iDESIGN) + '.csv', sep="\t")
