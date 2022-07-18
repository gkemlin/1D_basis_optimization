# Naive approach to build basis from the solutions of the quantum harmonic
# oscillator. We use
# - Nb local basis on each atom (LB)

include("include_common.jl")

output_dir="naive_results"
(!isdir(output_dir)) && (mkdir(output_dir))

if !isfile("$(output_dir)/λ_ref.json")
    data = Dict{String, Any}()
    # PDE parameters
    Nb_list = [1,2,3,4]
    Ng = 2001       # FD grid size inside box
    box_size = 20   # Box larger than I

    # Declare all general variables
    x_range, δx = discretize_space(Ng, box_size)
    a_sample = LinRange(0.1, 7, 150)

    data["Nb_list"] = Nb_list
    data["Ng"] = Ng
    data["box_size"] = box_size
    data["x_range"] = x_range
    data["δx"] = δx
    data["a_sample"] = a_sample

    λs_LB = Dict{Int64, Dict{Float64, Vector{Float64}}}()
    us_LB = Dict{Int64, Dict{Float64, Matrix{Float64}}}()
    cond_S_LB = Dict{Int64, Dict{Float64, Float64}}()

    # reference solution
    λs, us = exact_solution_all_a(a_sample, Ng, (x_range, δx))

    # generate variational solutions for different basis sets
    for Nb in Nb_list
        # global basis

        λs_LB_d = Dict{Float64, Vector{Float64}}()
        us_LB_d = Dict{Float64, Matrix{Float64}}()
        cond_S_LB_d = Dict{Float64, Float64}()

        for a in a_sample
            LBa = local_hermite_basis(a, Nb, x_range, δx)

            # Compute mass and overlap matrices
            Ha = - (1/2) * (1/δx)^2 * Δ(Ng) + Diagonal(V.(x_range, a))
            Sa = Symmetric(LBa'LBa)
            Mass_LBa = Symmetric(LBa'*(Ha*LBa))

            # variational solution in LBa
            λ_LB, C_LB = eigen(Mass_LBa, Sa)
            u_LB = LBa * C_LB[:,1:2]

            # update lists and dict
            λs_LB_d[a] = λ_LB[1:2]
            us_LB_d[a] = u_LB[:,1:2]
            cond_S_LB_d[a] = cond(Sa)
        end
        λs_LB[Nb] = λs_LB_d
        us_LB[Nb] = us_LB_d
        cond_S_LB[Nb] = cond_S_LB_d
    end
    open(io -> JSON3.write(io, data,  allow_inf=true), "naive_results/data.json", "w")
    open(io -> JSON3.write(io, λs,    allow_inf=true), "naive_results/λ_ref.json", "w")
    open(io -> JSON3.write(io, us,    allow_inf=true), "naive_results/u_ref.json", "w")
    open(io -> JSON3.write(io, λs_LB, allow_inf=true), "naive_results/λ_LB.json", "w")
    open(io -> JSON3.write(io, us_LB, allow_inf=true), "naive_results/u_LB.json", "w")
    open(io -> JSON3.write(io, cond_S_LB, allow_inf=true), "naive_results/cond_S_LB.json", "w")
end

# run plots
cd("$(output_dir)/")
include("$(output_dir)/plot_results.jl")
cd("../")
