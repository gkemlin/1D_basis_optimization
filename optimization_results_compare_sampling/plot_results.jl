# Plot results for the optimisation with the L criterion

include("../plots_common.jl")

refdir = pwd()

g_energy  = GroupPlot(3, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")
g_density = GroupPlot(3, 1, groupStyle="horizontal sep = 2cm, vertical sep = 1.75cm")

crit_list = [L"H^1{\rm-OBS}", L"L^2{\rm-OBS}", L"{\rm E-OBS}"]

for (i, dir) in enumerate(filter(d -> isdir(d), readdir()))

    display(dir)
    cd(joinpath(refdir, dir))

    data = open(JSON3.read, "data.json")
    Nb_list = data["Nb"]
    Ng = data["Ng"]
    box_size = data["box_size"]
    x_range = data["x_range"]
    δx = data["δx"]
    a_sample_list = data["a_sample_list"]
    a_sample_ext = Float64.(data["a_sample_ext"])

    λs    = open(JSON3.read, "λ_ref.json")
    λs_OB = open(JSON3.read, "λ_OB.json")
    λs_HB = open(JSON3.read, "λ_HB.json")
    us    = open(JSON3.read, "u_ref.json")
    us_OB = open(JSON3.read, "u_OB.json")
    us_HB = open(JSON3.read, "u_HB.json")

    cond_S_OB = open(JSON3.read, "cond_S_OB.json")
    cond_S_HB = open(JSON3.read, "cond_S_HB.json")

    # Plot
    # dissociation curve
    E = [sum(λs[a]) for a in a_sample_ext]
    E_OB = [[sum(λs_OB[ia_sample][a]) for a in a_sample_ext]
            for (ia_sample, a_sample) in enumerate(a_sample_list)]
    E_HB = [sum(λs_HB[1][a]) for a in a_sample_ext]

    p = [Plots.Linear(a_sample_ext, E_HB .-E, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, myorange"),
         Plots.Linear(a_sample_ext, E_OB[1] .- E, legendentry=latexstring("a \\in \\{1.5,5\\}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample_ext, E_OB[2] .- E, legendentry=latexstring(" a \\in \\{1.5,2.5\\}"),
                      mark="none", style="thick, dash dot, mygreen"),
         Plots.Linear(a_sample_ext, E_OB[3] .- E, legendentry=latexstring("a \\in \\{2.5\\}"),
                      mark="none", style="thick, dash dot dot, myviolet")]
    push!(g_energy, Axis(p, title=crit_list[i]*L"\quad {\rm Error\ on\ the\ energy}",
                         xlabel=L"$a$",
                         ylabel=L"\left|E_a - \widetilde{E}_a\right|", ymode="log",
                         legendStyle="at={(0.95,0.05)}, anchor=south east"))

    # Plot density errors
    ρs = []
    ρs_HB = []
    for a in a_sample_ext
        us_a    = reshape(us[a], (Ng,2)) ./ √δx
        us_HB_a = reshape(us_HB[1][a], (Ng,2)) ./ √δx
        push!(ρs,    us_a[:,1].^2 + us_a[:,2].^2)
        push!(ρs_HB, us_HB_a[:,1].^2 + us_HB_a[:,2].^2)
    end

    ρs_OB = []
    for (ia_sample, a_sample) in enumerate(a_sample_list)
        ρs_OB_sample = []
        for a in a_sample_ext
            us_OB_a = reshape(us_OB[ia_sample][a], (Ng,2)) ./ √δx
            push!(ρs_OB_sample, us_OB_a[:,1].^2 + us_OB_a[:,2].^2)
        end
        push!(ρs_OB, ρs_OB_sample)
    end

    err_HB = sum.([δx .* abs.(ρ - ρs_HB[i]) for (i,ρ) in enumerate(ρs)])
    err_OB = [sum.([δx .* abs.(ρ - ρs_OB[ia][i]) for (i,ρ) in enumerate(ρs)]) for ia in 1:3]

    p = [Plots.Linear(a_sample_ext, err_HB, legendentry=latexstring("\\text{HBS}"),
                      mark="none", style="thick, myorange"),
         Plots.Linear(a_sample_ext, err_OB[1], legendentry=latexstring("a \\in\\{1.5,5\\}"),
                      mark="none", style="thick, dashed, myblue"),
         Plots.Linear(a_sample_ext, err_OB[2], legendentry=latexstring("a \\in\\{1.5,2.5\\}"),
                      mark="none", style="thick, dash dot, mygreen"),
         Plots.Linear(a_sample_ext, err_OB[3], legendentry=latexstring("a \\in \\{2.5\\}"),
                      mark="none", style="thick, dash dot dot, myviolet")]
    push!(g_density, Axis(p, title=crit_list[i]*L"\quad {\rm Error\ on\ the\ density}",
                          xlabel=L"$x$",
                          ymode="log", ylabel=L"\Vert\rho_a - \widetilde\rho_a\Vert_{L^1}",
                          legendStyle="at={(0.95,0.05)}, anchor=south east"))
end
cd(refdir)
save("energy_compare_samples.pdf",  g_energy)
save("density_compare_samples.pdf", g_density)
