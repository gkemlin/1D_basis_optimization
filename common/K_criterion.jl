function K(R, offline_data, A)
    # Extract offline_data
    Ba_dict, Ma_dict = offline_data
    a_sample = Float64.(keys(Ba_dict))
    δa = compute_integral_weight(a_sample)

    # Compute integrand
    K_R = zero(Float64)
    I_R = IR(R)
    for a in a_sample
        Ba = Ba_dict[a]; Ma = Ma_dict[a]
        Sa = (Ba*I_R)'A*(Ba*I_R)
        K_R += tr(Ma*I_R*inv(Symmetric(Sa))*I_R')
    end
    # Multiply by Riemannian integral weight
    K_R * δa
end

"""
   Ba_dict and Ma_dict are the results of offline comutations
   (See offline_computations.jl)
"""
function ∇K(R, offline_data, A)
   # Precomputations
   Ba_dict, Ma_dict = offline_data
   N, Nb = size(R);
   a_sample = Float64.(keys(Ba_dict))
   δa = compute_integral_weight(a_sample)

   ∇K_loc_R = zero(R)
   I_R = IR(R)

   # loop over all a's
   for a in a_sample
       Ba = Ba_dict[a]; Ma = Ma_dict[a]
       # Precomputations
       Sa = (Ba*I_R)'A*(Ba*I_R)
       S_inv_a = inv(Symmetric(Sa))
       SB = Ba'A*Ba

       S1 = S_inv_a[1:Nb, 1:Nb]; S2 = S_inv_a[1:Nb, Nb+1:end]
       S3 = S_inv_a[Nb+1:end, Nb+1:end]

       SB1 = SB[1:N, 1:N]; SB2 = SB[1:N, N+1:end]
       SB3 = SB[N+1:end, N+1:end]

       Ma1 = Ma[1:N, 1:N]; Ma2 = Ma[1:N, N+1:end]
       Ma3 = Ma[N+1:end, N+1:end]



       ∇K_loc_R .+= 2 * (Ma1*R*S1 + Ma2*R*S2' + Ma2'R*S2 + Ma3*R*S3)
       ∇K_loc_R .-= 2 * (SB1 *R*S1*R' + SB2*R*S2'*R') * (Ma1 *R*S1 + Ma2*R*S2')
       ∇K_loc_R .-= 2 * (SB1 *R*S2*R' + SB2*R*S3 *R') * (Ma2'*R*S1 + Ma3*R*S2')
       ∇K_loc_R .-= 2 * (SB2'*R*S1*R' + SB3*R*S2'*R') * (Ma1 *R*S2 + Ma2*R*S3)
       ∇K_loc_R .-= 2 * (SB2'*R*S2*R' + SB3*R*S3 *R') * (Ma2'*R*S2 + Ma3*R*S3)
   end
   # Don't forget the Riemann integral weight !
   ∇K_loc_R *= δa
end
