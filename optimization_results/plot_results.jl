# Plot results for the optimisation with the L criterion

include("../plots_common.jl")
include("../common/finite_diff_solver.jl")

data = open(JSON3.read, "data.json")
Nb_list = data["Nb_list"]
Ng = data["Ng"]
N = data["N"]
box_size = data["box_size"]
x_range = data["x_range"]
δx = data["δx"]
a_sample = Float64.(data["a_sample_ext"])

ite_KL2 = open(JSON3.read, "ite_KL2.json")
ite_KH1 = open(JSON3.read, "ite_KH1.json")
ite_L   = open(JSON3.read, "ite_L.json")

OB_KL2 = open(JSON3.read, "OB_KL2.json")
OB_KH1 = open(JSON3.read, "OB_KH1.json")
OB_L   = open(JSON3.read, "OB_L.json")
HB     = open(JSON3.read, "HB.json")

λs        = open(JSON3.read, "λ_ref.json")
λs_OB_KL2 = open(JSON3.read, "λ_OB_KL2.json")
λs_OB_KH1 = open(JSON3.read, "λ_OB_KH1.json")
λs_OB_L   = open(JSON3.read, "λ_OB_L.json")
λs_HB     = open(JSON3.read, "λ_HB.json")

us        = open(JSON3.read, "u_ref.json")
us_OB_KL2 = open(JSON3.read, "u_OB_KL2.json")
us_OB_KH1 = open(JSON3.read, "u_OB_KH1.json")
us_OB_L   = open(JSON3.read, "u_OB_L.json")
us_HB     = open(JSON3.read, "u_HB.json")

crit_OB_KL2 = open(JSON3.read, "crit_OB_KL2.json")
crit_OB_KH1 = open(JSON3.read, "crit_OB_KH1.json")
crit_OB_L   = open(JSON3.read, "crit_OB_L.json")
crit_HB     = open(JSON3.read, "crit_HB.json")

cond_S_OB_KL2 = open(JSON3.read, "cond_S_OB_KL2.json")
cond_S_OB_KH1 = open(JSON3.read, "cond_S_OB_KH1.json")
cond_S_OB_L   = open(JSON3.read, "cond_S_OB_L.json")
cond_S_HB     = open(JSON3.read, "cond_S_HB.json")

# conditioning
g_cond_HB = Axis(title="Hermite basis",
                 xlabel=L"$a$",
                 ylabel="conditioning", ymode="log",
                 legendStyle="at={(0.95,0.95)}, anchor=north east")
g_cond_OB_KL2 = Axis(title="Criterion 1",
                     xlabel=L"$a$",
                     ylabel="conditioning", ymode="log",
                     legendStyle="at={(0.95,0.95)}, anchor=north east")
g_cond_OB_KH1 = Axis(title="Criterion 1",
                     xlabel=L"$a$",
                     ylabel="conditioning", ymode="log",
                     legendStyle="at={(0.95,0.95)}, anchor=north east")
g_cond_OB_L = Axis(title="Criterion 2",
                   xlabel=L"$a$",
                   ylabel="conditioning", ymode="log",
                   legendStyle="at={(0.95,0.95)}, anchor=north east")
style_list = ["thick, myred",
              "thick, dashed, myblue",
              "thick, dotted, mygreen",
              "thick, dash dot, myorange"]
for (i,Nb) in enumerate(Nb_list)
    style = style_list[i]
    cond_OB_KL2 = [cond_S_OB_KL2[Nb][a] for a in a_sample]
    cond_OB_KH1 = [cond_S_OB_KH1[Nb][a] for a in a_sample]
    cond_OB_L   = [cond_S_OB_L[Nb][a] for a in a_sample]
    cond_HB     = [cond_S_HB[Nb][a] for a in a_sample]
    push!(g_cond_HB, Plots.Linear(a_sample, cond_HB, legendentry=latexstring("N_b=$(Nb)"),
                                  mark="none", style=style))
    push!(g_cond_OB_KL2, Plots.Linear(a_sample, cond_OB_KL2, legendentry=latexstring("N_b=$(Nb)"),
                                      mark="none", style=style))
    push!(g_cond_OB_KH1, Plots.Linear(a_sample, cond_OB_KH1, legendentry=latexstring("N_b=$(Nb)"),
                                      mark="none", style=style))
    push!(g_cond_OB_L, Plots.Linear(a_sample, cond_OB_L, legendentry=latexstring("N_b=$(Nb)"),
                                    mark="none", style=style))
end
save("cond_S_HB.pdf", g_cond_HB)
save("cond_S_OB_KL2.pdf", g_cond_OB_KL2)
save("cond_S_OB_KH1.pdf", g_cond_OB_KH1)
save("cond_S_OB_L.pdf", g_cond_OB_L)

# process results
for Nb in Nb_list

    # dissociation curve
    E = [sum(λs[a]) for a in a_sample]
    E_OB_KL2 = [sum(λs_OB_KL2[Nb][a]) for a in a_sample]
    E_OB_KH1 = [sum(λs_OB_KH1[Nb][a]) for a in a_sample]
    E_OB_L   = [sum(λs_OB_L[Nb][a]) for a in a_sample]
    E_HB     = [sum(λs_HB[Nb][a]) for a in a_sample]

    g = GroupPlot(1, 2, groupStyle="vertical sep = 2cm")
    p = [Plots.Linear(a_sample, E, legendentry=L"{\rm reference}",
                      mark="none", style="thick, myred"),
         Plots.Linear(a_sample, E_OB_KL2, legendentry=latexstring("L^2\\text{-OBS}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, E_OB_KH1, legendentry=latexstring("H^1\\text{-OBS}"),
                      mark="none", style="thick, dash dot dot, myviolet"),
         Plots.Linear(a_sample, E_OB_L, legendentry=latexstring("\\text{E-OBS}"),
                      mark="none", style="thick, dotted, mygreen"),
         Plots.Linear(a_sample, E_HB, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, dash dot, myorange")]
    push!(g, Axis(p, title=latexstring("N_b=$(Nb)"),
                  xlabel=L"$a$", ylabel=L"\rm energy",
                  legendStyle="at={(0.95,0.05)}, anchor=south east"))
    p = [Plots.Linear(a_sample, E_OB_KL2.-E, legendentry=latexstring("L^2\\text{-OBS}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, E_OB_KH1.-E, legendentry=latexstring("H^1\\text{-OBS}"),
                      mark="none", style="thick, dash dot dot, myviolet"),
         Plots.Linear(a_sample, E_OB_L.-E, legendentry=latexstring("\\text{E-OBS}"),
                      mark="none", style="thick, dotted, mygreen"),
         Plots.Linear(a_sample, E_HB.-E, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, dash dot, myorange")]
    push!(g, Axis(p, title=latexstring("N_b=$(Nb)"),
                  xlabel=L"$a$",
                  ylabel=L"\left|E_a - \widetilde{E}_a\right|", ymode="log",
                  legendStyle="at={(0.95,0.95)}, anchor=north east"))
    save("Nb$(Nb)_energy.pdf", g)

    # densities
    ρs = []
    ρs_OB_KL2 = []
    ρs_OB_KH1 = []
    ρs_OB_L = []
    ρs_HB = []
    for a in a_sample
        us_a        = reshape(us[a], (Ng,2)) ./ √δx
        us_OB_KL2_a = reshape(us_OB_KL2[Nb][a], (Ng,2)) ./ √δx
        us_OB_KH1_a = reshape(us_OB_KH1[Nb][a], (Ng,2)) ./ √δx
        us_OB_L_a   = reshape(us_OB_L[Nb][a], (Ng,2)) ./ √δx
        us_HB_a     = reshape(us_HB[Nb][a], (Ng,2)) ./ √δx
        push!(ρs,    us_a[:,1].^2 + us_a[:,2].^2)
        push!(ρs_OB_KL2, us_OB_KL2_a[:,1].^2 + us_OB_KL2_a[:,2].^2)
        push!(ρs_OB_KH1, us_OB_KH1_a[:,1].^2 + us_OB_KH1_a[:,2].^2)
        push!(ρs_OB_L, us_OB_L_a[:,1].^2 + us_OB_L_a[:,2].^2)
        push!(ρs_HB, us_HB_a[:,1].^2 + us_HB_a[:,2].^2)
    end
    display([sum.(δx .* ρs) sum.(δx .* ρs_OB_KL2) sum.(δx .* ρs_OB_KH1) sum.(δx .* ρs_OB_L) sum.(δx .* ρs_HB)])

    id = Int(length(a_sample)/2)

    p = [Plots.Linear(x_range, ρs[id], legendentry=L"{\rm reference}",
                      mark="none", style="thick, myred"),
         Plots.Linear(x_range, ρs_OB_KL2[id], legendentry=latexstring("L^2\\text{-OBS}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(x_range, ρs_OB_KH1[id], legendentry=latexstring("H^1\\text{-OBS}"),
                      mark="none", style="thick, dash dot dot, myviolet"),
         Plots.Linear(x_range, ρs_OB_L[id], legendentry=latexstring("\\text{E-OBS}"),
                      mark="none", style="thick, dotted, mygreen"),
         Plots.Linear(x_range, ρs_HB[id], legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, dash dot, myorange")]
    g = Axis(p, title=latexstring("a=$(round(a_sample[Int(end/2)]; digits=2)), N_b=$(Nb)"),
             xlabel=L"$x$", xmin=-10, xmax=10,
             ylabel=L"\rho",
             legendStyle="at={(0.95,0.95)}, anchor=north east")
    save("Nb$(Nb)_density.pdf", g)

    # error L^1 norm
    g = GroupPlot(3, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")
    err_OB_KL2 = sum.([δx .* abs.(ρ - ρs_OB_KL2[i]) for (i,ρ) in enumerate(ρs)])
    err_OB_KH1 = sum.([δx .* abs.(ρ - ρs_OB_KH1[i]) for (i,ρ) in enumerate(ρs)])
    err_OB_L   = sum.([δx .* abs.(ρ - ρs_OB_L[i]) for (i,ρ) in enumerate(ρs)])
    err_HB     = sum.([δx .* abs.(ρ - ρs_HB[i]) for (i,ρ) in enumerate(ρs)])
    p = [Plots.Linear(a_sample, err_OB_KL2, legendentry=latexstring("L^2\\text{-OBS}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, err_OB_KH1, legendentry=latexstring("H^1\\text{-OBS}"),
                      mark="none", style="thick, dash dot dot, myviolet"),
         Plots.Linear(a_sample, err_OB_L, legendentry=latexstring("\\text{E-OBS}"),
                      mark="none", style="thick, dotted, mygreen"),
         Plots.Linear(a_sample, err_HB, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, dash dot, myorange")]
    push!(g, Axis(p, title=latexstring("N_b=$(Nb)"),
                  xlabel=L"$a$",
                  ymode="log", ylabel=L"\Vert\rho_a - \widetilde\rho_a\Vert_{L^1}",
                  legendStyle="at={(0.95,0.95)}, anchor=north east"))

    # error H^1 norm
    A_L2 = Diagonal(ones(Ng))
    A_H1 = A_L2 .- (1/δx^2)*Δ(Ng)
    err_OB_KL2 = sqrt.(sum.([δx .* (ρ - ρs_OB_KL2[i]) .* (A_H1 * (ρ - ρs_OB_KL2[i]))
                             for (i,ρ) in enumerate(ρs)]))
    err_OB_KH1 = sqrt.(sum.([δx .* (ρ - ρs_OB_KH1[i]) .* (A_H1 * (ρ - ρs_OB_KH1[i]))
                             for (i,ρ) in enumerate(ρs)]))
    err_OB_L   = sqrt.(sum.([δx .* (ρ - ρs_OB_L[i]) .* (A_H1 * (ρ - ρs_OB_L[i]))
                             for (i,ρ) in enumerate(ρs)]))
    err_HB     = sqrt.(sum.([δx .* (ρ - ρs_HB[i]) .* (A_H1 * (ρ - ρs_HB[i]))
                             for (i,ρ) in enumerate(ρs)]))
    p = [Plots.Linear(a_sample, err_OB_KL2, legendentry=latexstring("L^2\\text{-OBS}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, err_OB_KH1, legendentry=latexstring("H^1\\text{-OBS}"),
                      mark="none", style="thick, dash dot dot, myviolet"),
         Plots.Linear(a_sample, err_OB_L, legendentry=latexstring("\\text{E-OBS}"),
                      mark="none", style="thick, dotted, mygreen"),
         Plots.Linear(a_sample, err_HB, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, dash dot, myorange")]
    push!(g, Axis(p, title=latexstring("N_b=$(Nb)"),
                  xlabel=L"$a$",
                  ymode="log", ylabel=L"\Vert\rho_a - \widetilde\rho_a\Vert_{H^1}",
                  legendStyle="at={(0.95,0.95)}, anchor=north east"))

    # error von Weizsacker
    D = Tridiagonal(.-ones(Ng-1), zeros(Ng), ones(Ng-1)) / (2δx)
    err_OB_KL2 = sqrt.(abs.(sum.([δx .* (D*sqrt.(ρ) - D*sqrt.(ρs_OB_KL2[i])).^2
                                  for (i,ρ) in enumerate(ρs)])))
    err_OB_KH1 = sqrt.(abs.(sum.([δx .* (D*sqrt.(ρ) - D*sqrt.(ρs_OB_KH1[i])).^2
                                  for (i,ρ) in enumerate(ρs)])))
    err_OB_L   = sqrt.(abs.(sum.([δx .* (D*sqrt.(ρ) - D*sqrt.(ρs_OB_L[i])).^2
                                  for (i,ρ) in enumerate(ρs)])))
    err_HB     = sqrt.(abs.(sum.([δx .* (D*sqrt.(ρ) - D*sqrt.(ρs_HB[i])).^2
                                  for (i,ρ) in enumerate(ρs)])))
    p = [Plots.Linear(a_sample, err_OB_KL2, legendentry=latexstring("L^2\\text{-OBS}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, err_OB_KH1, legendentry=latexstring("H^1\\text{-OBS}"),
                      mark="none", style="thick, dash dot dot, myviolet"),
         Plots.Linear(a_sample, err_OB_L, legendentry=latexstring("\\text{E-OBS}"),
                      mark="none", style="thick, dotted, mygreen"),
         Plots.Linear(a_sample, err_HB, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, dash dot, myorange")]
    push!(g, Axis(p, title=latexstring("N_b=$(Nb)"),
                  xlabel=L"$a$",
                  ymode="log", ylabel=L"\Vert\nabla\sqrt{\rho_a} - \nabla\sqrt{\widetilde\rho_a}\Vert_{L^2}",
                  legendStyle="at={(0.95,0.95)}, anchor=north east"))


    save("Nb$(Nb)_density_error.pdf", g)
end

# plot functions
for Nb in Nb_list
    OB_KL2_Nb = reshape(OB_KL2[Nb], (Ng,Nb)) ./ √δx
    OB_KH1_Nb = reshape(OB_KH1[Nb], (Ng,Nb)) ./ √δx
    HB_Nb     = reshape(HB[Nb],     (Ng,N))    ./ √δx
    OB_L_Nb   = reshape(OB_L[Nb],   (Ng,Nb)) ./ √δx

    display(δx*OB_KL2_Nb'OB_KL2_Nb)
    display(δx*OB_KH1_Nb'OB_KH1_Nb)
    display(δx*OB_L_Nb'OB_L_Nb)
    display(δx*HB_Nb'HB_Nb)
    for μ in 1:Nb
        p = [Plots.Linear(x_range, OB_KL2_Nb[:,μ], legendentry=latexstring("L^2\\text{-OBS}"),
                          mark="none", style="thick, dashed, myblue"),
             Plots.Linear(x_range, OB_KH1_Nb[:,μ], legendentry=latexstring("H^1\\text{-OBS}"),
                          mark="none", style="thick, dash dot dot, myviolet"),
             Plots.Linear(x_range, OB_L_Nb[:,μ], legendentry=latexstring("\\text{E-OBS}"),
                          mark="none", style="thick, dotted, mygreen"),
             Plots.Linear(x_range, HB_Nb[:,μ], legendentry=L"{\rm HBS}",
                          mark="none", style="thick, dash dot, orange")]
        g = Axis(p, title=latexstring("\\mu=$(μ)"),
                 xlabel=L"$x$", xmin=-10, xmax=10,
                 legendStyle="at={(0.05,0.95)}, anchor=north west")
        save("X_$(μ)_Nb$(Nb).pdf",g)
    end
end

# print values of criterion
for Nb in Nb_list
    println(Nb)
    crOB_KL2 = crit_OB_KL2[Nb]
    crOB_KH1 = crit_OB_KH1[Nb]
    crOB_L   = crit_OB_L[Nb]
    crHB     = crit_HB[Nb]

    crit_list = ["KL2", "KH1", "L", "HB"]
    df = DataFrame(basis=crit_list)
    for crit in sort([c for c in keys(crOB_KL2)])
        df[!,crit] = [crOB_KL2[crit], crOB_KH1[crit], crOB_L[crit], crHB[crit]]
    end
    display(df)
end

# print convergence table
df = DataFrame(basis=["KL2", "KH1", "L"])
for Nb in Nb_list
    df[!,"$Nb"] = [ite_KL2[Nb], ite_KH1[Nb], ite_L[Nb]]
end
display(df)
