# Compute optimzed basis for both criteria.
# We then compare optimal basis (OB) for each criterion and localized hermite
# basis (HB).

include("include_common.jl")

output_dir="optimization_results"
(!isdir(output_dir)) && (mkdir(output_dir))

if !isfile("$(output_dir)/λ_ref.json")
    # PDE parameters
    Nb_list = [1,2,3,4]
    N = 10;         # Number of L2 basis functions to compute our basis functions
    Ng = 2001;      # FD grid size inside box
    box_size = 20   # computations are performed in [-box_size, box_size]

    # "a-sampling" parameters
    α = 1.5; β = 5
    Ns = 10
    Ns_ext = 100

    # tolerance for optimization
    tol = 1e-7

    # Declare all general variables
    x_range, δx = discretize_space(Ng, box_size)
    # range of a for the optimization
    a_sample = LinRange(α, β, Ns)
    # range of a for the dissociation curve
    a_sample_ext = LinRange(α, β, Ns_ext)

    data = Dict{String, Any}()

    data["Nb_list"] = Nb_list
    data["Ng"] = Ng
    data["N"] = N
    data["box_size"] = box_size
    data["x_range"] = x_range
    data["δx"] = δx
    data["a_sample"] = a_sample
    data["a_sample_ext"] = a_sample_ext

    # compute reference solution for a larger range a_sample_ext
    λs, us = exact_solution_all_a(a_sample_ext, Ng, (x_range, δx))

    # criterion K, L^2 norm
    A_L2 = Diagonal(ones(Ng))
    Ba_KL2, Ma_KL2 = run_offline_computations_K(a_sample, N, Ng, (x_range, δx); A=A_L2)
    offline_data_KL2 = Ba_KL2, Ma_KL2

    # criterion K, H^1 norm
    A_H1 = A_L2 .- (1/δx^2)*Δ(Ng)
    Ba_KH1, Ma_KH1 = run_offline_computations_K(a_sample, N, Ng, (x_range, δx); A=A_H1)
    offline_data_KH1 = Ba_KH1, Ma_KH1

    function build_optim_basis_K(Nb, offline_data, A)

        @info "Optimizing basis for criterion K with $(2Nb) basis"
        f(R) = -K(R, offline_data, A)
        function ∇f!(g, R)
            g .= -∇K(R, offline_data, A)
        end

        Rstart = R_init(N, Nb)
        res = optimize(R -> f(R), (g, R) -> ∇f!(g, R),
                       Rstart, LBFGS(manifold=Stiefel()),
                       Optim.Options(show_trace=true,
                                     g_tol=tol, iterations=500));

        res.minimizer, res.iterations
    end

    # criterion L
    Ba_L, Ma_L, λ_L, _ = run_offline_computations_L(a_sample, N, Ng, (x_range, δx))
    offline_data_L = Ba_L, Ma_L, λ_L

    function build_optim_basis_L(Nb, offline_data)

        @info "Optimizing basis for criterion L with $(2Nb) basis"
        function f∇f!(f, g, R)
            LR, λR, CR = L(R, offline_data)
            (g != nothing) && (g .= ∇L(R, λR, CR, offline_data))
            (f!= nothing) && (return LR)

        end

        Rstart = R_init(N, Nb)
        res = optimize(Optim.only_fg!(f∇f!),
                       Rstart, LBFGS(manifold=Stiefel()),
                       Optim.Options(show_trace=true,
                                     g_tol=tol, iterations=500));

        res.minimizer, res.iterations
    end

    # postprocess and save results
    ite_KL2_dict = Dict{Int64, Int64}()
    ite_KH1_dict = Dict{Int64, Int64}()
    ite_L_dict   = Dict{Int64, Int64}()

    λs_OB_KL2 = Dict{Int64, Dict{Float64, Vector{Float64}}}()
    λs_OB_KH1 = Dict{Int64, Dict{Float64, Vector{Float64}}}()
    λs_OB_L   = Dict{Int64, Dict{Float64, Vector{Float64}}}()
    λs_HB     = Dict{Int64, Dict{Float64, Vector{Float64}}}()

    us_OB_KL2 = Dict{Int64, Dict{Float64, Matrix{Float64}}}()
    us_OB_KH1 = Dict{Int64, Dict{Float64, Matrix{Float64}}}()
    us_OB_L   = Dict{Int64, Dict{Float64, Matrix{Float64}}}()
    us_HB     = Dict{Int64, Dict{Float64, Matrix{Float64}}}()

    crit_OB_KL2 = Dict{Int64, Dict{String, Float64}}()
    crit_OB_KH1 = Dict{Int64, Dict{String, Float64}}()
    crit_OB_L   = Dict{Int64, Dict{String, Float64}}()
    crit_HB     = Dict{Int64, Dict{String, Float64}}()

    cond_S_OB_KL2 = Dict{Int64, Dict{Float64, Float64}}()
    cond_S_OB_KH1 = Dict{Int64, Dict{Float64, Float64}}()
    cond_S_OB_L   = Dict{Int64, Dict{Float64, Float64}}()
    cond_S_HB     = Dict{Int64, Dict{Float64, Float64}}()

    OB_KL2_dict = Dict{Int64, Matrix{Float64}}()
    OB_KH1_dict = Dict{Int64, Matrix{Float64}}()
    OB_L_dict   = Dict{Int64, Matrix{Float64}}()
    HB_dict     = Dict{Int64, Matrix{Float64}}()

    # generate variational solutions for different basis sets
    for Nb in Nb_list

        # Launch optimization
        R_KL2, ite_KL2 = build_optim_basis_K(Nb, offline_data_KL2, A_L2)
        R_KH1, ite_KH1 = build_optim_basis_K(Nb, offline_data_KH1, A_H1)
        R_L, ite_L     = build_optim_basis_L(Nb, offline_data_L)
        R_HB           = R_init(N, Nb)

        # save iterations
        ite_KL2_dict[Nb] = ite_KL2
        ite_KH1_dict[Nb] = ite_KH1
        ite_L_dict[Nb]   = ite_L

        # save criterion values
        crit_OB_KL2_d = Dict{String, Float64}()
        crit_OB_KH1_d = Dict{String, Float64}()
        crit_OB_L_d   = Dict{String, Float64}()
        crit_HB_d     = Dict{String, Float64}()

        # KL2
        crit_OB_KL2_d["KL2"] = K(R_KL2, offline_data_KL2, A_L2)
        crit_OB_KL2_d["KH1"] = K(R_KL2, offline_data_KH1, A_H1)
        crit_OB_KL2_d["L"]   = L(R_KL2, offline_data_L)[1]

        # KH1
        crit_OB_KH1_d["KL2"] = K(R_KH1, offline_data_KL2, A_L2)
        crit_OB_KH1_d["KH1"] = K(R_KH1, offline_data_KH1, A_H1)
        crit_OB_KH1_d["L"]   = L(R_KH1, offline_data_L)[1]

        # L
        crit_OB_L_d["KL2"] = K(R_L, offline_data_KL2, A_L2)
        crit_OB_L_d["KH1"] = K(R_L, offline_data_KH1, A_H1)
        crit_OB_L_d["L"]   = L(R_L, offline_data_L)[1]

        # HB
        crit_HB_d["KL2"] = K(R_HB, offline_data_KL2, A_L2)
        crit_HB_d["KH1"] = K(R_HB, offline_data_KH1, A_H1)
        crit_HB_d["L"]   = L(R_HB, offline_data_L)[1]

        # save results anq QOI
        λs_OB_KL2_d = Dict{Float64, Vector{Float64}}()
        λs_OB_KH1_d = Dict{Float64, Vector{Float64}}()
        λs_OB_L_d = Dict{Float64, Vector{Float64}}()
        λs_HB_d   = Dict{Float64, Vector{Float64}}()

        us_OB_KL2_d = Dict{Float64, Matrix{Float64}}()
        us_OB_KH1_d = Dict{Float64, Matrix{Float64}}()
        us_OB_L_d = Dict{Float64, Matrix{Float64}}()
        us_HB_d   = Dict{Float64, Matrix{Float64}}()

        cond_S_OB_KL2_d = Dict{Float64, Float64}()
        cond_S_OB_KH1_d = Dict{Float64, Float64}()
        cond_S_OB_L_d = Dict{Float64, Float64}()
        cond_S_HB_d   = Dict{Float64, Float64}()

        for a in a_sample_ext
            HB = local_hermite_basis(a, Nb, x_range, δx)

            # optimized basis
            big_HB = local_hermite_basis(a, N, x_range, δx)
            OB_KL2 = big_HB*IR(R_KL2)
            OB_KH1 = big_HB*IR(R_KH1)
            OB_L = big_HB*IR(R_L)
            @assert size(HB) == size(OB_KL2)
            @assert size(HB) == size(OB_KH1)
            @assert size(HB) == size(OB_L)

            # discrete Hamiltonian
            Ha = - (1/2) * (1/δx)^2 * Δ(Ng) + Diagonal(V.(x_range, a))

            # solve for optimized basis K L2
            S_OB_KL2 = Symmetric(OB_KL2'OB_KL2)
            Mass_OB_KL2 = Symmetric(OB_KL2'*(Ha*OB_KL2))
            λ_OB_KL2, C_OB_KL2 = eigen(Mass_OB_KL2, S_OB_KL2)
            u_OB_KL2 = OB_KL2 * C_OB_KL2[:,1:2]
            # update lists and dict
            λs_OB_KL2_d[a]     = λ_OB_KL2[1:2]
            us_OB_KL2_d[a]     = u_OB_KL2[:,1:2]
            cond_S_OB_KL2_d[a] = cond(S_OB_KL2)

            # solve for optimized basis K H1
            S_OB_KH1 = Symmetric(OB_KH1'OB_KH1)
            Mass_OB_KH1 = Symmetric(OB_KH1'*(Ha*OB_KH1))
            λ_OB_KH1, C_OB_KH1 = eigen(Mass_OB_KH1, S_OB_KH1)
            u_OB_KH1 = OB_KH1 * C_OB_KH1[:,1:2]
            # update lists and dict
            λs_OB_KH1_d[a]     = λ_OB_KH1[1:2]
            us_OB_KH1_d[a]     = u_OB_KH1[:,1:2]
            cond_S_OB_KH1_d[a] = cond(S_OB_KH1)

            # solve for optimized basis L
            S_OB_L = Symmetric(OB_L'OB_L)
            Mass_OB_L = Symmetric(OB_L'*(Ha*OB_L))
            λ_OB_L, C_OB_L = eigen(Mass_OB_L, S_OB_L)
            u_OB_L = OB_L * C_OB_L[:,1:2]
            # update lists and dict
            λs_OB_L_d[a]     = λ_OB_L[1:2]
            us_OB_L_d[a]     = u_OB_L[:,1:2]
            cond_S_OB_L_d[a] = cond(S_OB_L)

            # solve for local hermite basis
            S_HB = Symmetric(HB'HB)
            Mass_HB = Symmetric(HB'*(Ha*HB))
            λ_HB, C_HB = eigen(Mass_HB, S_HB)
            u_HB = HB * C_HB[:,1:2]
            # update lists and dict
            λs_HB_d[a]     = λ_HB[1:2]
            us_HB_d[a]     = u_HB[:,1:2]
            cond_S_HB_d[a] = cond(S_HB)
        end
        λs_OB_KL2[Nb]     = λs_OB_KL2_d
        us_OB_KL2[Nb]     = us_OB_KL2_d
        crit_OB_KL2[Nb]   = crit_OB_KL2_d
        cond_S_OB_KL2[Nb] = cond_S_OB_KL2_d

        λs_OB_KH1[Nb]     = λs_OB_KH1_d
        us_OB_KH1[Nb]     = us_OB_KH1_d
        crit_OB_KH1[Nb]   = crit_OB_KH1_d
        cond_S_OB_KH1[Nb] = cond_S_OB_KH1_d

        λs_OB_L[Nb]     = λs_OB_L_d
        us_OB_L[Nb]     = us_OB_L_d
        crit_OB_L[Nb]   = crit_OB_L_d
        cond_S_OB_L[Nb] = cond_S_OB_L_d

        λs_HB[Nb]       = λs_HB_d
        us_HB[Nb]       = us_HB_d
        crit_HB[Nb]     = crit_HB_d
        cond_S_HB[Nb]   = cond_S_HB_d

        # centered basis functions
        big_HB = global_hermite_basis(N, x_range, δx)
        HB_dict[Nb]     = big_HB
        OB_KL2_dict[Nb] = big_HB*R_KL2
        OB_KH1_dict[Nb] = big_HB*R_KH1
        OB_L_dict[Nb]   = big_HB*R_L
    end
    # save results
    open(io -> JSON3.write(io, data,          allow_inf=true), "$(output_dir)/data.json", "w")
    open(io -> JSON3.write(io, λs,            allow_inf=true), "$(output_dir)/λ_ref.json", "w")
    open(io -> JSON3.write(io, us,            allow_inf=true), "$(output_dir)/u_ref.json", "w")

    open(io -> JSON3.write(io, λs_OB_KL2,     allow_inf=true), "$(output_dir)/λ_OB_KL2.json", "w")
    open(io -> JSON3.write(io, us_OB_KL2,     allow_inf=true), "$(output_dir)/u_OB_KL2.json", "w")
    open(io -> JSON3.write(io, OB_KL2_dict,   allow_inf=true), "$(output_dir)/OB_KL2.json", "w")
    open(io -> JSON3.write(io, ite_KL2_dict,  allow_inf=true), "$(output_dir)/ite_KL2.json", "w")
    open(io -> JSON3.write(io, crit_OB_KL2,   allow_inf=true), "$(output_dir)/crit_OB_KL2.json", "w")
    open(io -> JSON3.write(io, cond_S_OB_KL2, allow_inf=true), "$(output_dir)/cond_S_OB_KL2.json", "w")

    open(io -> JSON3.write(io, λs_OB_KH1,     allow_inf=true), "$(output_dir)/λ_OB_KH1.json", "w")
    open(io -> JSON3.write(io, us_OB_KH1,     allow_inf=true), "$(output_dir)/u_OB_KH1.json", "w")
    open(io -> JSON3.write(io, OB_KH1_dict,   allow_inf=true), "$(output_dir)/OB_KH1.json", "w")
    open(io -> JSON3.write(io, ite_KH1_dict,  allow_inf=true), "$(output_dir)/ite_KH1.json", "w")
    open(io -> JSON3.write(io, crit_OB_KH1,   allow_inf=true), "$(output_dir)/crit_OB_KH1.json", "w")
    open(io -> JSON3.write(io, cond_S_OB_KH1, allow_inf=true), "$(output_dir)/cond_S_OB_KH1.json", "w")

    open(io -> JSON3.write(io, λs_OB_L,       allow_inf=true), "$(output_dir)/λ_OB_L.json", "w")
    open(io -> JSON3.write(io, us_OB_L,       allow_inf=true), "$(output_dir)/u_OB_L.json", "w")
    open(io -> JSON3.write(io, OB_L_dict,     allow_inf=true), "$(output_dir)/OB_L.json", "w")
    open(io -> JSON3.write(io, ite_L_dict,    allow_inf=true), "$(output_dir)/ite_L.json", "w")
    open(io -> JSON3.write(io, crit_OB_L,     allow_inf=true), "$(output_dir)/crit_OB_L.json", "w")
    open(io -> JSON3.write(io, cond_S_OB_L,   allow_inf=true), "$(output_dir)/cond_S_OB_L.json", "w")

    open(io -> JSON3.write(io, λs_HB,         allow_inf=true), "$(output_dir)/λ_HB.json", "w")
    open(io -> JSON3.write(io, us_HB,         allow_inf=true), "$(output_dir)/u_HB.json", "w")
    open(io -> JSON3.write(io, HB_dict,       allow_inf=true), "$(output_dir)/HB.json", "w")
    open(io -> JSON3.write(io, crit_HB,       allow_inf=true), "$(output_dir)/crit_HB.json", "w")
    open(io -> JSON3.write(io, cond_S_HB,     allow_inf=true), "$(output_dir)/cond_S_HB.json", "w")
end

# run plots
cd("$(output_dir)/")
include("$(output_dir)/plot_results.jl")
cd("../")
