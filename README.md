Julia code for the optimization of basis in quantum chemistry on a toy model.
This code is used in [1] -> [insert arXiv link]


# Requirements:
Julia 1.7 with the libraries:
- Optim, Polynomials, SpecialPolynomials, KrylovKit for the computations;
- JSON3 for saving results;
- PGFPlots, LaTeXStrings, DataFrames for plotting results.

# Usage

This code runs the computations to perform the numerical experiments from [1].
Each section with numerical results is associated with a `.jl` file.
Run all the computations at once with `include("run_all_computations.jl")`.
You can also run computations for each section individually by
`include("file.jl")` to perform the simulations for the associated section.
It automatically saves results in the associated folder, in `.json` files.
Each folder contains a file `plot_results.jl` which generates `.pdf` files from
the data in the folder.

- `naive_approach.jl` runs simulations from Section 2.
- `optimize_basis.jl` runs simulations from Section 4.1.
- `optimize_basis_random_start.jl` runs simulations from Section 4.2.1.
- `optimize_basis_extrapolate_a.jl` runs simulations from Section 4.2.2.
- `optimize_basis_compare_sampling.jl` runs simulations from Section 4.2.3.
- `optimize_basis_N5.jl` runs simulations from Section 4.2.4.

# Contact
This is research code, not necessarily user-friendly or extremely robust.
If you have questions or encounter problems, get in touch!

