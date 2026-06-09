## Script that contains functions to fit trio-gcta models, either in the "classical" way or the AM-adjusted version, on simulated data from GNAMES and functions to fit all simulated datasets of one design in one go
# Sofie Perizonius

function GREMLModels.transform!(δ::Vector, θ::Vector, rmp, mode::Symbol)
    # Standard unconstrained trio GCTA transform
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[4]^2
    δ[3] = θ[3]^2 + θ[5]^2 + θ[6]^2
    δ[4] = θ[2] * θ[1]
    δ[5] = θ[3] * θ[1]
    δ[6] = θ[3] * θ[2] + θ[5] * θ[4]
    δ[7] = θ[7]
end

function fit_trio_gcta(gpath::String, mode::Symbol = :standard, rmp::Real = 0.0,
    Sigma::AbstractMatrix{<:Real} = I, residual::Real = 1.0, prints::Bool = false)

    # --- load data ---
    trio = SnpData(gpath)
    A = 2 * grm(trio.snparray)

    # --- index the individuals ---
    K = div(trio.people, 3)
    mid = 1:K
    pid = K+1:2K
    oid = 2K+1:3K

    # --- relationship matrices ---
    Amm = A[mid, mid]
    App = A[pid, pid]
    Aoo = A[oid, oid]
    Dpm = A[pid, mid] + A[mid, pid]
    Dom = A[oid, mid] + A[mid, oid]
    Dop = A[oid, pid] + A[pid, oid] 

    R = Diagonal(ones(K))
    r = [Amm, App, Aoo, Dpm, Dom, Dop, R]

    # --- fixed effects ---
    X = ones(K)

    # --- phenotype ---
    y = parse.(Float64, trio.person_info[oid, :phenotype])

    # --- data structure for VC model ---
    dat = GREMLData(y, X, r)
    
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

function save_model(model::GREMLModel, filename::String)
    JLD.save(filename, "model", model)
    println("Saved model to $filename")
end

function read_designs(csvfile::String)
    df = CSV.read(csvfile, DataFrame; header=false)

    # flatten everything into a vector of strings
    return vec(Matrix(df))
end

function run_design(iDESIGN:: Int64, designnames::Vector{<:AbstractString}, basedir=@__DIR__, outdir="results", 
    mode::Symbol = :standard, h2::Real = 0.0, Sigma::AbstractMatrix{<:Real} = I,
    residual::Real = 1.0, prints::Bool = false)
    println("Running design: $iDESIGN for $(length(designnames)) runs")
    println("Mode is $mode")
    results = DataFrame()
    fitting_results = DataFrame()
    thetas = DataFrame()
    se_fixed = DataFrame()
    se_vc = DataFrame()
    timings = DataFrame()

    rownames = ["Maternal indirect effect", "Paternal indirect effect", 
                "Offspring direct effect", "Maternal-paternal covariance",
               "Maternal-offspring covariance", "Paternal-offspring covariance",
                "Residual variance"]
                
    rownames_fitting = ["loglikelihood", "AIC", "BIC", "deviance"]
    insertcols!(results, 1, :rowname => rownames)
    insertcols!(fitting_results, 1, :rowname => rownames_fitting)

    
    # variance component names (same as your results)
    vc_names = rownames
    insertcols!(se_vc, 1, :rowname => vc_names)

    csvpath_results = joinpath(outdir, "Results_var_components_$(mode)_design_$(iDESIGN).csv")
    csvpath_timings = joinpath(outdir, "Timings_$(mode)_design_$(iDESIGN).csv")
    csvpath_fitting_results = joinpath(outdir, "Results_fitting_$(mode)_design_$(iDESIGN).csv")
    csvpath_thetas = joinpath(outdir, "Results_thetas_$(mode)_design_$(iDESIGN).csv")
    csvpath_se_fixed = joinpath(outdir, "Results_se_fixed_$(mode)_design_$(iDESIGN).csv")
    csvpath_se_vc = joinpath(outdir, "Results_se_vc_$(mode)_design_$(iDESIGN).csv")

    for designname in designnames
        println("\n==============================")
        println("Running design: $designname")
        println("==============================")

        # --- construct path to the SNP data file ---
        gpath = joinpath(
            basedir,
            "merged_reordered.$(designname)"
        )

        # --- check if file exists ---
        bedpath = "$(gpath).bed"

        if !isfile(bedpath)
            @warn "Skipping run $run — file not found: $gpath"
            continue
        end

        println("\n--- Fitting model for SETTING $designname ---")

        phenotypic_cor = match(r"RhoAM\.([0-9]+(?:\.[0-9]+)*)", designname)

        if phenotypic_cor !== nothing
            phenotypic_cor = phenotypic_cor.captures[1]  # This gives "0.00" as a string
            phenotypic_cor = parse(Float64, phenotypic_cor)
        end

        # genetic correlation between partners is heritability x phenotypic correlation
        genetic_cor = phenotypic_cor * h2

        #--- fit model ---
        model, t, _, _, _ = @timed fit_trio_gcta(gpath, mode, genetic_cor, Sigma, residual, prints)
        println(t)
        columnname = "$(designname)"

        timings[!, Symbol(columnname)] = [t]
        results[!, Symbol(columnname)] = model.δ
        fitting_results[!, Symbol(columnname)] = [loglikelihood(model), aic(model), bic(model), deviance(model)]
        thetas[!, Symbol(columnname)] = model.θ

        # fixed effects SE
        se_fixed[!, Symbol(columnname)] = GREMLModels.stderror(model)

        # variance component SE
        se_vc[!, Symbol(columnname)] = GREMLModels.stderrorvc(model)
        
        CSV.write(csvpath_results, results)
        CSV.write(csvpath_fitting_results, fitting_results)
        CSV.write(csvpath_thetas, thetas)
        CSV.write(csvpath_se_fixed, se_fixed)
        CSV.write(csvpath_se_vc, se_vc)
        CSV.write(csvpath_timings, timings)
    
        println("Completed RUN ✓")
    end

    println("\nAll runs completed.")
end
