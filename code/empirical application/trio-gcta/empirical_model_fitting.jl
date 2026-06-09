## Script used to fit the classical Trio-GCTA and AM-Trio-GCTA models to the empirical data
## Sofie Perizonius - 9-6-2026

using Pkg
Pkg.activate("Internship_2/vc_env")

using LinearAlgebra
using GREMLModels
using SnpArrays
using JLD
using StatsModels, StatsBase, Statistics
using DataFrames, CSV
using Serialization

include("empirical_trio_functions.jl")

# reload grms
bw_GRM = deserialize("empirical_trio/for_trio_gcta/bw_GRM.jls")
bl_GRM = deserialize("empirical_trio/for_trio_gcta/bl_GRM.jls")
height_GRM = deserialize("empirical_trio/for_trio_gcta/height_GRM.jls")
edu_GRM = deserialize("empirical_trio/for_trio_gcta/edu_GRM.jls")

# birth weight
sigma_bw = [
    0.1  0.00  0.00;
    0.00  0.00  0.00;
    0.00  0.00  0.3
    ]
residual_bw = 0.6

# standard
birth_weight_model_standard = fit_trio_gcta("empirical_trio/for_trio_gcta/birth_weight_trios", bw_GRM, sigma_bw, residual_bw, :standard, 0.00, false)
println("Results BW standard:")
show(birth_weight_model_standard)
print_full_δ(birth_weight_model_standard.δ)

# variance component SEse_vc[!, Symbol(columnname)] = GREMLModels.stderrorvc(model)
JLD.save("empirical_trio/for_trio_gcta/models/birth_weight_standard.jld", "model", birth_weight_model_standard)

#am_version
bw_phenotypic_cor = 0.037
bw_h2 = 0.3
bw_rmp = bw_phenotypic_cor * bw_h2

birth_weight_model_am = fit_trio_gcta("empirical_trio/for_trio_gcta/birth_weight_trios", bw_GRM, sigma_bw, residual_bw, :am_version, bw_rmp, false)

println("Results BW am:")
show(birth_weight_model_am)
print_full_δ(birth_weight_model_am.δ)

JLD.save("empirical_trio/for_trio_gcta/models/birth_weight_am.jld", "model", birth_weight_model_am)

# birth length
sigma_bl = [
    0.1  0.00  0.00;
    0.00  0.00  0.00;
    0.00  0.00  0.15
    ]
residual_bl = 0.8

# standard
birth_length_model_standard = fit_trio_gcta("empirical_trio/for_trio_gcta/birth_length_trios", bl_GRM, sigma_bl, residual_bl, :standard, 0.00, false)
println("Results BL standard:")
show(birth_length_model_standard)
print_full_δ(birth_length_model_standard.δ)

JLD.save("empirical_trio/for_trio_gcta/models/birth_length_standard.jld", "model", birth_length_model_standard)

#am_version
bl_phenotypic_cor = 0.154
bl_h2 = 0.13
bl_rmp = bl_phenotypic_cor * bl_h2

birth_length_model_am = fit_trio_gcta("empirical_trio/for_trio_gcta/birth_length_trios", bl_GRM, sigma_bl, residual_bl, :am_version, bl_rmp, false)
println("Results BL am:")
show(birth_length_model_am)

print_full_δ(birth_length_model_am.δ)

JLD.save("empirical_trio/for_trio_gcta/models/birth_length_am.jld", "model", birth_length_model_am)

# height
sigma_height = [
    0.05  0.00  0.00;
    0.00  0.05  0.00;
    0.00  0.00  0.3
    ]
residual_height = 0.6

# standard
height_model_standard = fit_trio_gcta("empirical_trio/for_trio_gcta/height_trios", height_GRM, sigma_height, residual_height, :standard, 0.00, false)
println("Results height standard:")
show(height_model_standard)
print_full_δ(height_model_standard.δ)

JLD.save("empirical_trio/for_trio_gcta/models/height_standard.jld", "model", height_model_standard)

#am_version
height_phenotypic_cor = 0.256
height_h2 = 0.55
height_rmp = height_phenotypic_cor * height_h2

height_model_am = fit_trio_gcta("empirical_trio/for_trio_gcta/height_trios", height_GRM, sigma_height, residual_height, :am_version, height_rmp, false)

println("Results height am:")
show(height_model_am)
print_full_δ(height_model_am.δ)

JLD.save("empirical_trio/for_trio_gcta/models/height_am.jld", "model", height_model_am)

# cito scores
sigma_edu = [
    0.05  0.00  0.00;
    0.00  0.05  0.00;
    0.00  0.00  0.15
    ]
residual_edu = 0.8

# standard
edu_model_standard = fit_trio_gcta("empirical_trio/for_trio_gcta/edu_trios", edu_GRM, sigma_edu, residual_edu, :standard, 0.00, false)

println("Results edu standard:")
show(edu_model_standard)
print_full_δ(edu_model_standard.δ)

JLD.save("empirical_trio/for_trio_gcta/models/edu_standard.jld", "model", edu_model_standard)

#am_version
edu_phenotypic_cor = 0.426
edu_h2 = 0.20
edu_rmp = edu_phenotypic_cor * edu_h2

edu_model_am = fit_trio_gcta("empirical_trio/for_trio_gcta/edu_trios", edu_GRM, sigma_edu, residual_edu, :am_version, edu_rmp, false)

println("Results edu am:")
show(edu_model_am)
print_full_δ(edu_model_am.δ)

JLD.save("empirical_trio/for_trio_gcta/models/edu_am.jld", "model", edu_model_am)

## standard errors
birth_length_model_standard = JLD.load("empirical_trio/for_trio_gcta/models/birth_length_standard.jld")["model"]
birth_length_model_am = JLD.load("empirical_trio/for_trio_gcta/models/birth_length_am.jld")["model"]

hessian!(birth_length_model_standard.opt.H, birth_length_model_standard)
V_bl_standard = inv(0.5 * birth_length_model_standard.opt.H)
sqrt.(max.(LinearAlgebra.diag(V_bl_standard), 0.00))

stderrorvc(birth_length_model_standard)

hessian!(birth_length_model_am.opt.H, birth_length_model_am)
stderrorvc(birth_length_model_am)

birth_weight_model_standard = JLD.load("empirical_trio/for_trio_gcta/models/birth_weight_standard.jld")["model"]
birth_weight_model_am = JLD.load("empirical_trio/for_trio_gcta/models/birth_weight_am.jld")["model"]

hessian!(birth_weight_model_standard.opt.H, birth_weight_model_standard)

V_bw_standard = inv(0.5 * birth_weight_model_standard.opt.H)
sqrt.(max.(LinearAlgebra.diag(V_bw_standard), 0.00))

stderrorvc(birth_weight_model_standard)

hessian!(birth_weight_model_am.opt.H, birth_weight_model_am)
stderrorvc(birth_weight_model_am)

height_model_standard = JLD.load("empirical_trio/for_trio_gcta/models/height_standard.jld")["model"]
height_model_am = JLD.load("empirical_trio/for_trio_gcta/models/height_am.jld")["model"]

hessian!(height_model_standard.opt.H, height_model_standard)
stderrorvc(height_model_standard)

hessian!(height_model_am.opt.H, height_model_am)
stderrorvc(height_model_am)

edu_model_standard = JLD.load("empirical_trio/for_trio_gcta/models/edu_standard.jld")["model"]
edu_model_am = JLD.load("empirical_trio/for_trio_gcta/models/edu_am.jld")["model"]

hessian!(edu_model_standard.opt.H, edu_model_standard)
V_edu_standard = inv(0.5 * edu_model_standard.opt.H)
sqrt.(max.(LinearAlgebra.diag(V_edu_standard), 0.00))

stderrorvc(edu_model_standard)

hessian!(edu_model_am.opt.H, edu_model_am)
stderrorvc(edu_model_am)
