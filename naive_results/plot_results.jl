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
λs_GB = open(JSON3.read, "λ_GB.json")
λs_LB = open(JSON3.read, "λ_LB.json")
us    = open(JSON3.read, "u_ref.json")
us_GB = open(JSON3.read, "u_GB.json")
us_LB = open(JSON3.read, "u_LB.json")

cond_S_LB = open(JSON3.read, "cond_S_LB.json")

g_cond = Axis(xlabel=L"$a$",
              ylabel=L"\rm conditioning", ymode="log",
              legendStyle="at={(0.95,0.95)}, anchor=north east")
style_list = ["thick, dashed, myblue",
              "thick, dotted, mygreen"]
for (i,Nb) in enumerate(Nb_list)
    style = style_list[i]
    # conditioning
    cond_S = [cond_S_LB[Nb][a] for a in a_sample]
    push!(g_cond, Plots.Linear(a_sample, cond_S, legendentry=latexstring("N_b=$(2Nb)"),
                               mark="none", style=style))
end
cond_1 = [1.0 for a in a_sample]
push!(g_cond, Plots.Linear(a_sample, cond_1, legendentry=latexstring("1"),
                           mark="none", style="thick, black"))
save("cond_S_LB.pdf", g_cond)

# process results
for Nb in Nb_list

    # dissociation curve
    E = [sum(λs[a]) for a in a_sample]
    E_GB = [sum(λs_GB[Nb][a]) for a in a_sample]
    E_LB = [sum(λs_LB[Nb][a]) for a in a_sample]

    g = GroupPlot(2, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")
    p = [Plots.Linear(a_sample, E, legendentry=L"{\rm reference}",
                      mark="none", style="thick, myred"),
         Plots.Linear(a_sample, E_GB, legendentry=latexstring("{\\rm global}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, E_LB, legendentry=latexstring("{\\rm local}"),
                      mark="none", style="thick, dotted, mygreen")]
    push!(g, Axis(p, title=latexstring("N_b=$(2Nb)"),
                  xlabel=L"$a$", ylabel=L"\rm energy", ymin=0.5, ymax=3,
                  legendStyle="at={(0.05,0.95)}, anchor=north west"))
    p = [Plots.Linear(a_sample, E_GB.-E, legendentry=latexstring("{\\rm global}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, E_LB.-E, legendentry=latexstring("{\\rm local}"),
                      mark="none", style="thick, dotted,  mygreen")]
    push!(g, Axis(p, title=latexstring("N_b=$(2Nb)"),
                  xlabel=L"$a$",
                  ylabel=L"\left|E_a - \widetilde{E}_a\right|", ymode="log",
                  legendStyle="at={(0.05,0.95)}, anchor=north west"))
    save("Nb$(2Nb)_energy.pdf", g)

    # densities
    ρs = []
    ρs_GB = []
    ρs_LB = []
    for a in a_sample
        us_a    = reshape(us[a], (Ng,2)) ./ √δx
        us_GB_a = reshape(us_GB[Nb][a], (Ng,2)) ./ √δx
        us_LB_a = reshape(us_LB[Nb][a], (Ng,2)) ./ √δx
        push!(ρs,    us_a[:,1].^2 + us_a[:,2].^2)
        push!(ρs_GB, us_GB_a[:,1].^2 + us_GB_a[:,2].^2)
        push!(ρs_LB, us_LB_a[:,1].^2 + us_LB_a[:,2].^2)
    end
    display([sum.(δx .* ρs) sum.(δx .* ρs_GB) sum.(δx .* ρs_LB)])
    err_GB = sum.([δx .* abs.(ρ - ρs_GB[i]) for (i,ρ) in enumerate(ρs)])
    err_LB = sum.([δx .* abs.(ρ - ρs_LB[i]) for (i,ρ) in enumerate(ρs)])

    id = Int(length(a_sample)/2)
    g = GroupPlot(2, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")
    p = [Plots.Linear(x_range, ρs[id], legendentry=L"{\rm reference}",
                      mark="none", style="thick, myred"),
         Plots.Linear(x_range, ρs_GB[id], legendentry=latexstring("{\\rm global}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(x_range, ρs_LB[id], legendentry=latexstring("{\\rm local}"),
                      mark="none", style="thick, dotted, mygreen")]
    push!(g, Axis(p, title=latexstring("a=$(round(a_sample[Int(end/2)]; digits=2)), N_b=$(2Nb)"),
                  xlabel=L"$x$", xmin=-10, xmax=10,
                  ylabel=L"\rho",
                  legendStyle="at={(0.05,0.95)}, anchor=north west"))
    p = [Plots.Linear(a_sample, err_GB, legendentry=latexstring("{\\rm global}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample, err_LB, legendentry=latexstring("{\\rm local}"),
                      mark="none", style="thick, dotted, mygreen")]
    push!(g, Axis(p, title=latexstring("N_b=$(2Nb)"),
                  xlabel=L"$a$",
                  ymode="log", ylabel=L"\Vert\rho_a - \widetilde\rho_a\Vert_{L^1}",
                  legendStyle="at={(0.95,0.05)}, anchor=south east"))
    save("Nb$(2Nb)_density.pdf", g)
end
