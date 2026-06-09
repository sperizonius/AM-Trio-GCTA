## Script that contains functions to fit trio-gcta models, either in the "classical" way or the AM-adjusted version, 
## on simulated data from GNAMES and functions to fit all simulated datasets of one design in one go

## Difference with functions used on simulated data is that these functions that a precomputed grm and take the (residualized) phenotype information 
## from a csv file

## Sofie Perizonius - 9-6-2026

function GREMLModels.transform!(δ::Vector, θ::Vector, rmp, mode::Symbol)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[4]^2
    δ[3] = θ[3]^2 + θ[5]^2 + θ[6]^2
    δ[4] = θ[2] * θ[1]
    δ[5] = θ[3] * θ[1]
    δ[6] = θ[3] * θ[2] + θ[5] * θ[4]
    δ[7] = θ[7]
end

function fit_trio_gcta(gpath::String, grm::Matrix{Float64}, Sigma::AbstractMatrix{<:Real} = I, residual::Real = 1.0, mode::Symbol = :standard, rmp::Real = 0.0, prints::Bool = true)
    # --- load data ---
    trio = SnpData(gpath)

    # --- index the individuals ---
    K = div(trio.people, 3)
    mid = 1:K
    pid = K+1:2K
    oid = 2K+1:3K

    A = grm

    # --- relationship matrices ---
    Amm = A[mid, mid]
    App = A[pid, pid]
    Aoo = A[oid, oid]
    Dpm = A[pid, mid] + A[mid, pid]
    Dom = A[oid, mid] + A[mid, oid]
    Dop = A[oid, pid] + A[pid, oid] 

    R = Diagonal(ones(K))
    r = [Amm, App, Aoo, Dpm, Dom, Dop, R]

    # --- covariates ---
    covariates_phenotypes = CSV.read(gpath * ".csv", DataFrame; header=true)

    X = ones(K)
    y = covariates_phenotypes[oid, 4]
    y = parse.(Float64, y)

    println("Variance in y is:")
    println(var(y))

    println("Mean is")
    println(mean(y))

    # --- data structure for VC model ---
    dat = GREMLData(y, X, r)

    # --- initial δ and bounds ---
    epsilon = 1e-6
    Sigma_pd = Sigma + epsilon * I  # I is the identity operator

    L = cholesky(Sigma_pd).L  # extract lower triangular factor
    
    t1 = L[1,1]
    t2 = L[2,1]
    t3 = L[3,1]
    t4 = L[2,2]
    t5 = L[3,2]
    t6 = L[3,3]

    # --- initial δ and bounds ---
    δ_init  = [t1, t2, t3, t4, t5, t6, residual]
    δ_lower = [0.0, -Inf, -Inf, 0.0, -Inf, 0.0, 0.0]

    # --- fit model ---
    println("Fitting the model....")
    m = GREMLModel(dat, δ_init, δ_lower; mode = mode, rmp = rmp)
    
    @time fit!(m; verbose=prints)

    return m
end

function print_full_δ(δ::Vector)
    labels = [
        "Maternal indirect effect",
        "Paternal indirect effect",
        "Offspring direct effect",
        "Maternal-paternal covariance",
        "Maternal-offspring covariance",
        "Paternal-offspring covariance",
        "Residual variance"
    ]
    for i in 1:length(δ)
        println("$(labels[i]) = $(δ[i])")
    end
end

# function to save model after fitting
function save_model(model::GREMLModel, filename::String)
    JLD.save(filename, "model", model)
    println("Saved model to $filename")
end

# function for loading saved models
function load_model(filename::String)
    data = JLD.load(filename)
    return data["model"]
end
