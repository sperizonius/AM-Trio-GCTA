## Script to make reorder files that only contain the IDs of the children to be used with PLINK to make PLINK files for the children
## Sofie Perizonius - 9-6-2026

# make children subset files
bl_reorder = pd.read_csv("for_trio_gcta/reorder_birth_length.txt", sep=r"\s+", header=None)
n_bl = int(bl_reorder.shape[0] /3)
bl_reorder_children = bl_reorder.iloc[2*n_bl : 3*n_bl + 1, :]
bl_reorder_children.to_csv('for_trio_gcta/birth_length_children.txt', sep="\t", header=False, index=False)

bw_reorder = pd.read_csv("for_trio_gcta/reorder_birth_weight.txt", sep=r"\s+", header=None)
n_bw = int(bw_reorder.shape[0] /3)
bw_reorder_children = bw_reorder.iloc[2*n_bw : 3*n_bw + 1, :]
bw_reorder_children.to_csv('for_trio_gcta/birth_weight_children.txt', sep="\t", header=False, index=False)

height_reorder = pd.read_csv("for_trio_gcta/reorder_height.txt", sep=r"\s+", header=None)
n_height = int(height_reorder.shape[0] /3)
height_reorder_children = height_reorder.iloc[2*n_height : 3*n_height + 1, :]
height_reorder_children.to_csv('for_trio_gcta/height_children.txt', sep="\t", header=False, index=False)

edu_reorder = pd.read_csv("for_trio_gcta/reorder_edu.txt", sep=r"\s+", header=None)
n_edu = int(edu_reorder.shape[0] /3)
edu_reorder_children = edu_reorder.iloc[2*n_edu : 3*n_edu + 1, :]
edu_reorder_children.to_csv('for_trio_gcta/edu_children.txt', sep="\t", header=False, index=False)
