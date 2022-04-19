@doc raw"""
    Compute new criterion ``L_I(R) = \int_I |λ_FD(a) - λ(a,R)|^2 da``
    λ_FD is a dictionnary of eigenvalues computed with finite difference
    over the sample a_sample.
"""
function L(R, offline_data)
    Ba_dict, Ma_dict, λ_FD = offline_data
    a_sample = Float64.(keys(Ba_dict))

    # Compute all λ(a,R)
    λ = Dict(); C = Dict()
    I_R = IR(R)
    for a in a_sample
        # Compute overlap
        Ba = Ba_dict[a]
        S = (Ba*I_R)'*(Ba*I_R)
        inv_sqrt_S = sqrt(inv(Symmetric(S)))
        # Diagonalize Hc
        Hc_aR = I_R'Ma_dict[a]*I_R
        D, U, info = eigsolve(Symmetric(inv_sqrt_S'*Hc_aR*inv_sqrt_S),
                              2, :SR)
        λ[a] = D[1:2]; C[a] = inv_sqrt_S*[U[1] U[2]];
        # Test conditioning
        (cond(S)>1e6) && (@warn "cond S" cond(S))
        Pa = C[a]*C[a]';
        (norm(Pa*S*Pa .- Pa)>1e-10) && (@warn "test Pa" norm(Pa*S*Pa .- Pa))
    end

    # Compute integral
    δa = compute_integral_weight(a_sample)
    L_R = sum( abs2(sum(λ_FD[a]) - sum(λ[a])) for a in a_sample ) * δa
    return L_R, λ, C
end

function ∇E_R(R, a, Ma_dict, C)
    Ma = Ma_dict[a]; Ca_tilda = C[a];
    Pa = Ca_tilda*Ca_tilda'
    n1 = div(size(Ma,1),2); n2 = div(size(Pa,1),2)
    ∇E = Ma[1:n1,1:n1]*R*Pa[1:n2,1:n2] .+
        Ma[1:n1,n1+1:end]*R*Pa[n2+1:end,1:n2] .+
        Ma[n1+1:end,1:n1]*R*Pa[1:n2,n2+1:end] .+
        Ma[n1+1:end,n1+1:end]*R*Pa[n2+1:end,n2+1:end]
    2 * ∇E
end

function ∇E_C(R, a, Ba_dict, λ, C)
    Λa = diagm(λ[a]); Ca = C[a]; Ba = Ba_dict[a]
    Qa = Ca*Λa*Ca'
    Sa = Ba'Ba
    # Sum all blocs
    n1 = div(size(Sa,1),2); n2 = div(size(Qa,1),2)
    ∇E = Sa[1:n1,1:n1]*R*Qa[1:n2,1:n2] .+
        Sa[1:n1,n1+1:end]*R*Qa[n2+1:end,1:n2] .+
        Sa[n1+1:end,1:n1]*R*Qa[1:n2,n2+1:end] .+
        Sa[n1+1:end,n1+1:end]*R*Qa[n2+1:end,n2+1:end]
    - 2 * ∇E
end

∇E(R, a, Ba_dict, Ma_dict, λ_FD, λ, C) = ∇E_R(R, a, Ma_dict, C) .+
    ∇E_C(R, a, Ba_dict, λ, C)

function ∇L(R, λ, C,  offline_data)
    Ba_dict, Ma_dict, λ_FD = offline_data
    a_sample = Float64.(keys(Ba_dict))

    # First compute the integrand
    #  Intruction without ForwardDiff. To be used after debug.
    ∇L_R = sum((sum(λ_FD[a]) - sum(λ[a])) .* ∇E(R, a, Ba_dict, Ma_dict, λ_FD, λ, C) for
                a in a_sample)

    # Compute weight
    δa = compute_integral_weight(a_sample)
    ∇L_R .* (-2*δa)
end
