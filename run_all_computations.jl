# Run all the simulations
# -----------------------

# Run simulations from Section 2.
include("naive_approach.jl")

# Run simulations from Section 4.1.
include("optimize_basis.jl")

# Run simulations from Section 4.2.1.
include("optimize_basis_random_start.jl")

# Run simulations from Section 4.2.2.
include("optimize_basis_extrapolate_a.jl")

# Run simulations from Section 4.2.3.
include("optimize_basis_compare_sampling.jl")

# Run simulations from Section 4.2.4.
include("optimize_basis_N5.jl")

