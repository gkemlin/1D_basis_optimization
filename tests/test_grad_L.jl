using PyPlot
include("../include_common.jl")

# PDE parameters
N = 10;         # Number of L2 basis functions
Ng = 1501;      # FD grid size inside box
box_size = 30   # Box larger than I

# Declare all general variables
x_range, δx = discretize_space(Ng, box_size)

a_sample = LinRange(1,5,15)
R = R_init(N,4)
H = rand(N,4); H /= norm(H);

Ba_dict, Ma_dict, λ_FD, u_FD = run_offline_computations_L(a_sample, N, Ng, (x_range, δx))
offline_data = (Ba_dict, Ma_dict, λ_FD)

# test grad
L_R, λ, C = L(R, offline_data)
∇L_R = ∇L(R, λ, C, offline_data)


list_norms = []
list_test = []
list_ε = LinRange(1e-5,1e-1,100)

y = 1
h = 1

for ε in list_ε
    L_RεH, _ = L(R .+ ε*H, offline_data)
    yεh = (y + ε*h)^2
    push!(list_norms, norm(L_RεH .- L_R .- ε*tr(∇L_R'H)))
    push!(list_test, norm(yεh - y^2 - 2*y*ε*h))
end

figure()
loglog(list_ε, list_norms, label="NOUS")
loglog(list_ε, list_test, label="x^2")
legend()
