using PyPlot
include("../include_common.jl")

# PDE parameters
N = 10;         # Number of L2 basis functions
Ng = 1501;      # FD grid size inside box
box_size = 30   # Box larger than I

# Declare all general variables
x_range, δx = discretize_space(Ng, box_size)

# norm for K
A_L2 = diagm(ones(Ng))
A_H1 = A_L2 .- (1/δx^2)*Tridiagonal(ones(Ng-1), -2*ones(Ng), ones(Ng-1))
A = A_H1

a_sample = LinRange(1,5,15)
R = R_init(N,4)
H = rand(N,4); H /= norm(H);

offline_data = run_offline_computations_K(a_sample, N, Ng, (x_range, δx); A)

# test for a given a
#  Ba_dict, Ma_dict = offline_data
#  Ba = Ba_dict[3.0]
#  Ma = Ma_dict[3.0]
#  Σa_B = Ba[:,1:N]'Ba[:,N+1:end]
#  I_R = IR(R)
#  S_R = (Ba*I_R)'*(Ba*I_R)
#  S_inv_R = inv(Symmetric(S_R))
#  K_R = tr(I_R*S_inv_R*I_R'*Ma)
#  Q_R = (I_R*S_inv_R)' * Ma * (I_R*S_inv_R)

#  list_norms = []
#  list_test = []
#  list_ε = LinRange(1e-5,1e-1,100)

#  y = 1
#  h = 1

#  for ε in list_ε
#      RεH = R .+ ε .* H
#      I_RεH = IR(RεH)
#      I_H = IR(H)
#      S_RεH = (Ba*I_RεH)'*(Ba*I_RεH)
#      S_inv_RεH = inv(Symmetric(S_RεH))
#      K_RεH = tr(I_RεH*S_inv_RεH*I_RεH'*Ma)
#      HΣR = H'Σa_B*R + R'Σa_B*H
#      mid = vcat( hcat(H'R + R'H, HΣR), hcat(HΣR', H'R + R'H) )
#      dK_R = 2*tr(I_H * S_inv_RεH * I_R' * Ma)
#      dK_R -= tr(Q_R*mid)
#      push!(list_norms, norm(K_RεH .- K_R .- ε * dK_R))
#      yεh = (y + ε*h)^2
#      push!(list_test, norm(yεh - y^2 - 2*y*ε*h))
#  end

#  figure()
#  loglog(list_ε, list_norms, label="NOUS")
#  loglog(list_ε, list_test, label="x^2")
#  legend()

# test grad
K_R = K(R, offline_data, A)
∇K_R = ∇K(R, offline_data, A)


list_norms = []
list_test = []
list_ε = LinRange(1e-5,1e-1,100)

y = 1
h = 1

for ε in list_ε
    K_RεH = K(R .+ ε*H, offline_data, A)
    yεh = (y + ε*h)^2
    push!(list_norms, norm(K_RεH .- K_R .- ε*tr(∇K_R'H)))
    push!(list_test, norm(yεh - y^2 - 2*y*ε*h))
end

figure()
loglog(list_ε, list_norms, label="NOUS")
loglog(list_ε, list_test, label="x^2")
legend()
