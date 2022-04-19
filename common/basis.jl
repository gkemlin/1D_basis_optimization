"""
Hermite-polynomial Gaussian type basis. Return as functions of x.
If parity is set to "even" (resp. "odd"), returns only the N first
even (resp. odd) functions.
"""
function simple_hermite_basis(N)
    X = variable(Polynomial{Rational{Int}})
    c(n) = 1.0 / √(2^(n-1) * factorial(big(n-1)) * √π)

    # All basis function
    hermite_pol = [basis(Hermite, n)(X) for n in 0:2*N]
    [x->p(x)*exp(-(x)^2/2) * c(n) for (n,p) in enumerate(hermite_pol)][1:N]
end

"""
Same as above but the basis functions are evaluated on the finite difference grid.
The result is an array of size Ng×N.
"""
function global_hermite_basis(N, x_range, δx)
    Ng = length(x_range)
    # Analytical basis
    hermite_basis = simple_hermite_basis(N)
    # Discretized basis
    B = zeros(Ng, N);
    for i in 1:N
        B[:,i] = hermite_basis[i].(x_range)
    end
    B .* √δx
end

"""
Translate the simple basis given by the above routine on -a and a.
The result is an Array of size Ng×2N.
"""
function local_hermite_basis(a, N, x_range, δx)
    DB_on_grid = zeros(length(x_range), 2*N)
    centered_basis = simple_hermite_basis(N)
    # Place basis element in - and + a
    for i in 1:N
        DB_on_grid[:,i] .= centered_basis[i].(x_range .- a)
        DB_on_grid[:,N+i] .= centered_basis[i].(x_range .+ a)
    end
    DB_on_grid .* √δx
end
