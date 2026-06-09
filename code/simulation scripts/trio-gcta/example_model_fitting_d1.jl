## Example script for fitting AM-Trio-GCTA models to all simulated datasets for design 1

using Pkg
Pkg.activate("../vc_env")

using LinearAlgebra
using GREMLModels
using SnpArrays
using JLD
using StatsModels, StatsBase, Statistics
using DataFrames, CSV

include("fitting_trio_functions.jl")

## DESIGN_1
design_1 = read_designs(joinpath(dirname(@__DIR__), "simulated_data/DESIGN_1/to_analyze/setting_list_DESIGN.1_am_update.csv"))
Sigma_1 = [
    0.00  0.00  0.00;
    0.00  0.00  0.00;
    0.00  0.00  0.6
    ]

run_design(1, design_1, joinpath(dirname(@__DIR__), "simulated_data/DESIGN_1/to_analyze"), 
            joinpath(dirname(@__DIR__), "simulation_results/DESIGN_1/results_am_version"), 
            :am_version, 0.6, Sigma_1, 0.4, false)
