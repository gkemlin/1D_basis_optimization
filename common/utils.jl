"""
Function concerning creation and handling of R matrix
"""
function R_init(N, Nb; random_start=false)
    (Nb ≤ 0) && (error("Must set Nb > 0"))
    if random_start
        qr(Float64.(randn(N, Nb))).Q[:,1:Nb]
    else
        Float64.(Matrix(I, N, Nb))
    end
end

function rand_unitary_matrix(N)
    M = rand(N,N)
    exp((M .- M') ./ 2)
end
R_init_rand(N, Nb) = rand_unitary_matrix(N)*R_init(N,Nb)

IR(R) = vcat( hcat(R, zero(R)), hcat(zero(R), R))


"""
Compute integral weigth in K and L criterion.
For now only returns the weight for a regular sampling.
"""
function compute_integral_weight(a_sample)
    a_sample_sorted = sort(a_sample)
    (length(a_sample)≠1) && (return a_sample_sorted[2] - a_sample_sorted[1])
    one(Float64)
end

