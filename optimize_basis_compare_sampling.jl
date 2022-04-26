# Compute optimzed basis for both criterion
# We thenc compare optimal basis (OB) and localized hermite basis (HB)

include("include_common.jl")

# PDE parameters
Nb = 3
N = 10;         # Number of L2 basis functions to compute our basis functions
Ng = 2001;      # FD grid size inside box
box_size = 20   # Box larger than I
x_range, δx = discretize_space(Ng, box_size)
tol = 1e-7

# Criterion optimization
A_L2 = Diagonal(ones(Ng))
A_H1 = A_L2 .- (1/δx^2)*Δ(Ng)
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

    res.minimizer
end

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

    res.minimizer
end

# Output dir
ref_dir = "optimization_results_compare_sampling"
for (i, dir) in enumerate(["Crit1_L2", "Crit1_H1", "Crit2"])
    output_dir = joinpath(ref_dir, dir)
    (!isdir(output_dir)) && (mkdir(output_dir))

    if !isfile("$(output_dir)/λ_ref.json")
        # "a-sampling" parameters
        α = 1.5; β = 5
        Ns_ext = 100

        # range of a for the optimization
        a_sample_list = [[1.5,5], [1.5, 2.5], [2.5]]
        # range of a for the dissociation curve
        a_sample_ext = LinRange(α, β, Ns_ext)

        data = Dict{String, Any}()

        data["Nb"] = Nb
        data["Ng"] = Ng
        data["N"] = N
        data["box_size"] = box_size
        data["x_range"] = x_range
        data["δx"] = δx
        data["a_sample_list"] = a_sample_list
        data["a_sample_ext"] = a_sample_ext

        # compute reference solution for a larger range a_sample_ext
        λs, us = exact_solution_all_a(a_sample_ext, Ng, (x_range, δx))

        # postprocess and save results
        λs_OB = Dict{Int64, Dict{Float64, Vector{Float64}}}()
        us_OB = Dict{Int64, Dict{Float64, Matrix{Float64}}}()
        crit_OB = Dict{Int64, Dict{String, Float64}}()
        cond_S_OB = Dict{Int64, Dict{Float64, Float64}}()
        OB_dict = Dict{Int64, Matrix{Float64}}()

        λs_HB     = Dict{Int64, Dict{Float64, Vector{Float64}}}()
        us_HB     = Dict{Int64, Dict{Float64, Matrix{Float64}}}()
        crit_HB     = Dict{Int64, Dict{String, Float64}}()
        cond_S_HB     = Dict{Int64, Dict{Float64, Float64}}()
        HB_dict     = Dict{Int64, Matrix{Float64}}()

        # generate variational solutions for different basis sets
        for (ia_sample, a_sample) in enumerate(a_sample_list)
            if i == 1 # crit1 L^2
                A = A_L2
                offline_data = run_offline_computations_K(a_sample, N, Ng, (x_range, δx); A)
                R_opti = build_optim_basis_K(Nb, offline_data, A)
            elseif i == 2 # crit1 H^1
                A = A_H1
                offline_data = run_offline_computations_K(a_sample, N, Ng, (x_range, δx); A)
                R_opti = build_optim_basis_K(Nb, offline_data, A)
            else # crit2
                A = A_L2
                offline_data = run_offline_computations_L(a_sample, N, Ng, (x_range, δx))
                R_opti = build_optim_basis_L(Nb, offline_data)
            end

            # save results anq QOI
            λs_OB_d = Dict{Float64, Vector{Float64}}()
            λs_HB_d   = Dict{Float64, Vector{Float64}}()

            us_OB_d = Dict{Float64, Matrix{Float64}}()
            us_HB_d   = Dict{Float64, Matrix{Float64}}()

            cond_S_OB_d = Dict{Float64, Float64}()
            cond_S_HB_d   = Dict{Float64, Float64}()

            for a in a_sample_ext
                HB = local_hermite_basis(a, Nb, x_range, δx)

                # optimized basis
                big_HB = local_hermite_basis(a, N, x_range, δx)
                OB = big_HB*IR(R_opti)
                @assert size(HB) == size(OB)

                # discrete Hamiltonian
                Ha = - (1/2) * (1/δx)^2 * Δ(Ng) + Diagonal(V.(x_range, a))

                # solve for optimized basis w.r. to chosen criteria
                S_OB = Symmetric(OB'OB)
                Mass_OB = Symmetric(OB'*(Ha*OB))
                λ_OB, C_OB = eigen(Mass_OB, S_OB)
                u_OB = OB * C_OB[:,1:2]

                # update lists and dict
                λs_OB_d[a]     = λ_OB[1:2]
                us_OB_d[a]     = u_OB[:,1:2]
                cond_S_OB_d[a] = cond(S_OB)

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
            λs_OB[ia_sample]     = λs_OB_d
            us_OB[ia_sample]     = us_OB_d
            cond_S_OB[ia_sample] = cond_S_OB_d

            λs_HB[ia_sample]       = λs_HB_d
            us_HB[ia_sample]       = us_HB_d
            cond_S_HB[ia_sample]   = cond_S_HB_d

            # centered basis functions
            big_HB = global_hermite_basis(N, x_range, δx)
            HB_dict[ia_sample]     = big_HB
            OB_dict[ia_sample]     = big_HB*R_opti
        end
        # save results
        open(io -> JSON3.write(io, data,          allow_inf=true), "$(output_dir)/data.json", "w")
        open(io -> JSON3.write(io, λs,            allow_inf=true), "$(output_dir)/λ_ref.json", "w")
        open(io -> JSON3.write(io, us,            allow_inf=true), "$(output_dir)/u_ref.json", "w")

        open(io -> JSON3.write(io, λs_OB,     allow_inf=true), "$(output_dir)/λ_OB.json", "w")
        open(io -> JSON3.write(io, us_OB,     allow_inf=true), "$(output_dir)/u_OB.json", "w")
        open(io -> JSON3.write(io, OB_dict,   allow_inf=true), "$(output_dir)/OB.json", "w")
        open(io -> JSON3.write(io, crit_OB,   allow_inf=true), "$(output_dir)/crit_OB.json", "w")
        open(io -> JSON3.write(io, cond_S_OB, allow_inf=true), "$(output_dir)/cond_S_OB.json", "w")

        open(io -> JSON3.write(io, λs_HB,         allow_inf=true), "$(output_dir)/λ_HB.json", "w")
        open(io -> JSON3.write(io, us_HB,         allow_inf=true), "$(output_dir)/u_HB.json", "w")
        open(io -> JSON3.write(io, HB_dict,       allow_inf=true), "$(output_dir)/HB.json", "w")
        open(io -> JSON3.write(io, crit_HB,       allow_inf=true), "$(output_dir)/crit_HB.json", "w")
        open(io -> JSON3.write(io, cond_S_HB,     allow_inf=true), "$(output_dir)/cond_S_HB.json", "w")
    end
end

# run plots
cd("$(ref_dir)/")
include("$(ref_dir)/plot_results.jl")
cd("../")
