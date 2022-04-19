# Toy model to solve
#
# -1/2Δu + Vu = λu
#
# with V(x) = (1/(8a^2+1)) (x-a)^2(x+a)^2

V(x, a) = (1/(8a^2 + 4)) * (x-a)^2 * (x+a)^2

# resolution by finite differences on a big domain (5*a size) with dirichlet boundary
# conditions Ng points
function discretize_space(Ng, box_size)
    x_range = LinRange(-box_size, box_size, Ng)
    δx = x_range[2] - x_range[1]
    (x_range, δx)
end

Δ(n) = Tridiagonal(ones(n-1), -2.0*ones(n), ones(n-1))
function exact_solution(a, Ng, FD_grid)
    println("-------------------------------------")
    println("Exact solution for a=$(a)")
    # space discretisation
    x_range, δx = FD_grid

    # discrete hamiltonian -1/2 Δ + V
    H = - (1/2) * (1/δx)^2 * Δ(Ng) + Diagonal(V.(x_range, a))

    # solve the eigenproblem
    vals, vecs, info = eigsolve(H, 2, :SR; issymmetric=true,
                                verbosity=1, tol=1e-12, maxiter=1000,
                                eager=true)
    @show(size(vecs[1]))
    println("Eigenvalues for a=$(a):")
    display(vals)

    (vals, hcat(vecs...))
end

"""
reference solution for a given range of a
"""
function exact_solution_all_a(a_sample, Ng, FD_grid)
    λa_dict = Dict(); Ψa_dict = Dict();
    for a in a_sample
        # Solve problem by finite difference
        λa, Ψa = exact_solution(a, Ng, FD_grid)
        λa_dict[a] = λa[1:2]
        Ψa_dict[a] = Ψa[:,1:2]
    end
    λa_dict, Ψa_dict
end

