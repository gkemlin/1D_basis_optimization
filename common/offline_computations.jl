"""
K
"""
function run_offline_computations_K(a_sample, N, Ng, FD_grid;
                                    A=diagm(ones((length(FD_grid[1])))),
                                    full_output=false)

    @info "Offline computations for criterion K"
    # Precomputations
    Ba_dict = Dict(); Ma_dict = Dict();
    λa_dict = Dict(); Ψa_dict = Dict();
    x_range, δx = FD_grid

    for a in a_sample
        # Solve problem by finite difference
        λa, Ψa = exact_solution(a, Ng, FD_grid)
        # Assemble dicts
        Ba = local_hermite_basis(a, N, x_range, δx)
        Pa = (Ψa[:,1]*Ψa[:,1]' + Ψa[:,2]*Ψa[:,2]')
        @assert norm(Pa*Pa .- Pa) < 1e-10
        Ba_dict[a] = Ba
        Ma_dict[a] = Ba'*A*Pa*A*Ba
        λa_dict[a] = λa[1:2]
        Ψa_dict[a] = Ψa[:,1:2]
    end
    if full_output
        Ba_dict, Ma_dict, λa_dict, Ψa_dict
    else
        Ba_dict, Ma_dict
    end
end

"""
L
"""
function run_offline_computations_L(a_sample, N, Ng, FD_grid)
    @info "Offline computations for criterion L"
    # Precomputations
    Ba_dict = Dict(); Ma_dict = Dict();
    λa_dict = Dict(); Ψa_dict = Dict();
    x_range, δx = FD_grid

    for a in a_sample
        # Solve problem by finite difference
        λa, Ψa = exact_solution(a, Ng, FD_grid)
        Ha = - (1/2) * (1/δx)^2 * Δ(Ng) + Diagonal(V.(x_range, a))
        Ba = local_hermite_basis(a, N, x_range, δx)
        Ma_dict[a] = Ba'Ha*Ba
        Ba_dict[a] = Ba
        λa_dict[a] = λa[1:2]
        Ψa_dict[a] = Ψa[:,1:2]
    end
    Ba_dict, Ma_dict, λa_dict, Ψa_dict
end
