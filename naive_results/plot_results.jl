# Plot results for the naive approach

include("../plots_common.jl")

data = open(JSON3.read, "data.json")
Nb_list = data["Nb_list"]
Ng = data["Ng"]
box_size = data["box_size"]
x_range = data["x_range"]
δx = data["δx"]
a_sample = data["a_sample"]

λs    = open(JSON3.read, "λ_ref.json")
λs_LB = open(JSON3.read, "λ_LB.json")
us    = open(JSON3.read, "u_ref.json")
us_LB = open(JSON3.read, "u_LB.json")

cond_S_LB = open(JSON3.read, "cond_S_LB.json")

g_cond = Axis(xlabel=L"$a$", xmode="log",
              ylabel=L"\rm conditioning", ymode="log",
              legendStyle="at={(0.95,0.95)}, anchor=north east")
style_list = ["thick, dashed, myblue", "thick, dotted, mygreen",
              "thick, dash dot, myorange", "thick, dash dot dot, myviolet"]
for (i,Nb) in enumerate(Nb_list)
    style = style_list[i]
    # conditioning
    cond_S = [cond_S_LB[Nb][a] for a in a_sample]
    push!(g_cond, Plots.Linear(a_sample, cond_S, legendentry=latexstring("N_b=$(Nb)"),
                               mark="none", style=style))
end
cond_1 = [1.0 for a in a_sample]
push!(g_cond, Plots.Linear(a_sample, cond_1, legendentry=latexstring("1"),
                           mark="none", style="thick, black"))
save("cond_S_LB.pdf", g_cond)

# process results
# Reference data
E = [sum(λs[a]) for a in a_sample]
ρs = []
for a in a_sample
    us_a    = reshape(us[a], (Ng,2)) ./ √δx
    push!(ρs,    us_a[:,1].^2 + us_a[:,2].^2)
end

# Data for hermite basis
E_LB = []
ρs_LB = []
err_LB = []

for Nb in Nb_list
    # dissociation curve
    push!(E_LB, [sum(λs_LB[Nb][a]) for a in a_sample])

    # densities
    ρs_LB_Nb = []
    for a in a_sample
        us_LB_a = reshape(us_LB[Nb][a], (Ng,2)) ./ √δx
        push!(ρs_LB_Nb, us_LB_a[:,1].^2 + us_LB_a[:,2].^2)
    end
    push!(err_LB, sum.([δx .* abs.(ρ - ρs_LB_Nb[i]) for (i,ρ) in enumerate(ρs)]))
    push!(ρs_LB, ρs_LB_Nb)
end

g = GroupPlot(2, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")

# Plot dissociation curves
p = [Plots.Linear(a_sample, E, legendentry=L"{\rm reference}",
                  mark="none", style="thick, myred")]
for (i, Nb) in enumerate(Nb_list)
    push!(p, Plots.Linear(a_sample, E_LB[i], legendentry=latexstring("{\\rm HBS:\\;}N_b=$(Nb)"),
                          mark="none", style=style_list[i]))
end
push!(g, Axis(p; title=latexstring(" "),
              xlabel=L"$a$", ylabel=L"\rm energy", ymin=0.5, ymax=2.5,
              legendStyle="at={(0.95,0.95)}, anchor=north east"))
# Plot errors on dissociation curves
p = [Plots.Linear(a_sample, E_LB[i] .- E,
                  legendentry=latexstring("{\\rm HBS:\\;}N_b=$(Nb)"),
                  mark="none", style=style_list[i])
     for (i, Nb) in enumerate(Nb_list)]
push!(g, Axis(p; title=latexstring(" "),
              xlabel=L"$a$",
              ylabel=L"\left|E_a - \widetilde{E}_a\right|", ymode="log",
              legendStyle="at={(0.95,0.95)}, anchor=north east"))
save("dissociation_curves_hermite.pdf", g)

# plot densities and errors on densities
id = Int(length(a_sample)/2)
g = GroupPlot(2, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")
p = [Plots.Linear(x_range, ρs[id], legendentry=L"{\rm reference}",
                  mark="none", style="thick, myred")]
for (i, Nb) in enumerate(Nb_list)
    push!(p, Plots.Linear(x_range, ρs_LB[i][id],
                 legendentry=latexstring("{\\rm HBS:\\;}N_b=$(Nb)"),
                 mark="none", style=style_list[i]))
end
push!(g, Axis(p, title=latexstring("a=$(round(a_sample[Int(end/2)]; digits=2))"),
              xlabel=L"$x$", xmin=-10, xmax=10,
              ylabel=L"\rho",
              legendStyle="at={(0.05,0.95)}, anchor=north west"))
p = [Plots.Linear(a_sample, err_LB[i],
                  legendentry=latexstring("{\\rm HBS:\\;}N_b=$(Nb)"),
                  mark="none", style=style_list[i])
     for (i, Nb) in enumerate(Nb_list)]
push!(g, Axis(p, title=" ",
              xlabel=L"$a$",
              ymode="log", ylabel=L"\Vert\rho_a - \widetilde\rho_a\Vert_{L^1}",
              legendStyle="at={(0.95,0.05)}, anchor=south east"))
save("density_hermite.pdf", g)
