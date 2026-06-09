## Example script for computation of the GRM as input for the Trio-GCTA models
## Sofie Perizonius - 9-6-2026

using Pkg
Pkg.activate("Internship_2/vc_env")

using LinearAlgebra
using GREMLModels
using SnpArrays
using JLD
using StatsModels, StatsBase, Statistics
using DataFrames, CSV
using JLD2
using Serialization

include("empirical_trio_functions.jl")

## Compute GRMs and save them
# height
height_gpath = "empirical_trio/for_trio_gcta/height_trios"
@time height_trio = SnpData(height_gpath)
@time height_A = 2 * grm(height_trio.snparray)

# save
serialize("empirical_trio/for_trio_gcta/height_GRM.jls", height_A)

