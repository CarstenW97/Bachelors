using ColorSchemes, CairoMakie, FileIO, JuMP
using ProgressMeter

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel

imgpath = joinpath("docs", "imgs", "EMP-model")

#=
Warm up problem to get all the symbols in the model.
Create arrays using these symbols as the names.
This is the fast way of doing:
mu = Float64[]
atp = Float64[]
etc.
=#
test_model = GlnModel.gln_model();
syms = Symbol.(all_variables(test_model))
for sym in syms
    eval(:($sym=Float64[]))
end

prev_sol=Dict{Symbol, Float64}() # a variable to store the previous solution

ratio_lb = 0.05
ratio_ub = 10
@showprogress for var in reverse(range(log(ratio_lb); stop=log(ratio_ub), length=20))
    nadh_nad_ratio = exp(var)

    gln_model = GlnModel.gln_model(
        glc_ext = 50e-3,
        lac_ext = 1e-4,
        nh4_ext = 10e-3,
        ac_ext = 1e-16,
        etoh_ext = 1e-4,
        # atp_adp_ratio = -1, # unconstrained
        nadh_nad_ratio = nadh_nad_ratio,
        num_ms = 10,
    )

    solve_1 = GlnModel.max_mu!(gln_model, prev_sol)

    if !solve_1
        continue # if mu not found, skip
    end

    solve_2 = GlnModel.find_nearest!(gln_model, prev_sol) # find a solution close to the previous one, smoothen things

    if solve_1 && solve_2
        for sym in syms
            push!(eval(sym), value(gln_model[sym]))
        end
        # update old previous solution
        GlnModel.update!(gln_model, prev_sol)
    end
end

function plot_proteome(x, fracs, fraclabels; legendlabel="",xlabel="", ylabel="", xscale=identity)
    # NB: can't display more than 10 and not less than 3
    z = cumsum(fracs'; dims = 1) # no zeros
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    N = length(fraclabels)
    cs = eval(Meta.parse("ColorSchemes.Spectral_$N"))

    f = Figure()
    ax = Axis(f[1, 1], xscale=xscale, xlabel=xlabel, ylabel=ylabel)

    for i = 1:N
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i], color=cs[i])
    end

    f[1, 2] = Legend(f, ax, legendlabel, unique = true, framevisible = false)

    return f
end

#=
Plot cofactor metabolites
=#
fracs = exp.(hcat([adp, atp, nad, nadh]...))
fraclabels = ["adp", "atp", "nad", "nadh"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Metabolites",
    xlabel="NADH/NAD ratio",
    ylabel="Relative abundance",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Cofactor_Metabolites_NADH_NAD.pdf"), f)

#=
Plot upper glycolysis metabolites
=#
fracs = exp.(hcat([g6p, f6p, fdp, dhap, g3p, dpg13]...))
fraclabels = ["g6p", "f6p", "fdp", "dhap", "g3p", "dpg13"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Metabolites",
    xlabel="NADH/NAD ratio",
    ylabel="Relative abundance",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Metabolites_NADH_NAD.pdf"), f)

#=
Plot lower glycolysis metabolites
=#
fracs = exp.(hcat([pg3, pg2, pep]...))
fraclabels = ["pg3", "pg2", "pep"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Metabolites",
    xlabel="NADH/NAD ratio",
    ylabel="Relative abundance",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Lower_Glycolysis_Metabolites_NADH_NAD.pdf"), f)

#=
Plot Gibbs dissipation of reactions in upper glycolysis
=#
fracs = hcat([dg_pts.*v_pts, dg_pgi .* v_pgi, dg_pfk .* v_pfk, dg_fba .*v_fba, dg_tpi .* v_tpi, dg_gapd .* v_gapd]...)
# fracs = hcat([dg_pts, dg_pgi, dg_pfk, dg_fba, dg_tpi, dg_gapd]...)
fraclabels = ["pts", "pgi", "pfk", "fba", "tpi", "gapd"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Reaction",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Gibbs dissipation",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Gibbs_Dissapation_NADH_NAD.pdf"), f)


#=
Plot upper glycolysis proteome
=#
fracs = hcat([pts, pgi, pfk, fba, tpi, gapd]...)
fraclabels = ["pts", "pgi", "pfk", "fba", "tpi", "gapd"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Enzyme",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Enzyme concentration",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Proteome_NADH_NAD.pdf"), f)

#=
Plot lower glycolysis proteome
=#
fracs = hcat([pgk, pgm, pyk, eno, ldh]...)
fraclabels = ["pgk", "pgm", "pyk", "eno", "ldh"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Enzyme",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Enzyme concentration",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Lower_Glycolysis_Proteome_NADH_NAD.pdf"), f)

#=
Plot glutamine synthesis proteome
=#
fracs = hcat([ppc, pfl, pdh, cs, aconta, acontb, icdh, gludy, glns]...)
fraclabels = ["ppc", "pfl", "pdh", "cs", "aconta", "acontb", "icdh", "gludy", "glns"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Enzyme",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Enzyme concentration",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Glutamine_synthesis_Proteome_NADH_NAD.pdf"), f)

#=
Plot remaining proteome (transporters; ethanol -> accoa; acetate -> accoa)
=#
fracs = hcat([lact, nh4t, fort, etoht, alcd, acald, act, ackr, ptar]...)
fraclabels = ["lact", "nh4t", "fort", "etoht", "alcd", "acald", "act", "ackr", "ptar"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Enzyme",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Enzyme concentration",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Remaining_Proteome_NADH_NAD.pdf"), f)

#=
Plot Fluxes of reactions in upper glycolysis
=#
fracs = hcat([v_pts,v_pgi, v_pfk, v_fba, v_tpi, v_gapd]...)
fraclabels = ["v_pts", "v_pgi", "v_pfk", "v_fba", "v_tpi", "v_gapd"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Reaction",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Flux",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Fluxes_NADH_NAD.pdf"), f)

#=
Plot Fluxes of reactions in glutamine synthesis
=#
fracs = hcat([v_ppc, v_pfl, v_pdh, v_cs, v_aconta, v_acontb, v_icdh, v_gludy, v_glns]...)
fraclabels = ["v_ppc", "v_pfl", "v_pdh", "v_cs", "v_aconta", "v_acontb", "v_icdh", "v_gludy", "v_glns"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Reaction",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Flux",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Glutamine_synthesis_Fluxes_NADH_NAD.pdf"), f)

#=
Plot dG of reactions in upper glycolysis
=#

fracs = hcat([dg_pts, dg_pgi, dg_pfk, dg_fba, dg_tpi, dg_gapd]...)
fraclabels = ["dG_pts", "dG_pgi", "dG_pfk", "dG_fba", "dG_tpi", "dG_gapd"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Reaction",
    xlabel="NADH/NAD ratio",
    ylabel="Relative dG of reaction",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_dG_NADH_NAD.pdf"), f)

#=
Plot dG of reactions in glycolysis synthesis
=#
fracs = hcat([dg_ppc, dg_pfl, dg_pdh, dg_cs, dg_aconta, dg_acontb, dg_icdh, dg_gludy, dg_glns]...)
fraclabels = ["dg_ppc", "dg_pfl", "dg_pdh", "dg_cs", "dg_aconta", "dg_aconta", "dg_icdh", "dg_gludy", "dg_glns"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Reaction",
    xlabel="NADH/NAD ratio",
    ylabel="Relative dG of reaction",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Glutamine_synthesis_dG_NADH_NAD.pdf"), f)

### Absolte plots ###

#=
Plot cofactor metabolites
=#
x = exp.(nadh)./exp.(nad)
a = exp.([adp]...)
b = exp.([atp]...)
c = exp.([nad]...)
d = exp.([nadh]...)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d)
line4 = lines!(ax, x, d)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "Metabolite concentration [M]"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="Metabolite",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4]],
    ["ADP", "ATP", "NAD", "NADH"],
)
f

FileIO.save(joinpath(imgpath, "Cofactor_Metabolites_NADH_NAD_abs.pdf"), f)

#=
Plot upper glycolysis metabolites
=#
x = exp.(nadh)./exp.(nad)
a = exp.([g6p]...)
b = exp.([f6p]...)
c = exp.([fdp]...)
d = exp.([dhap]...)
e = exp.([g3p]...)
g = exp.([dpg13]...)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d)
line4 = lines!(ax, x, d)
scatter5 = scatter!(ax, x, e)
line5 = lines!(ax, x, e)
scatter6 = scatter!(ax, x, g)
line6 = lines!(ax, x, g)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "Metabolite concentration [M]"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="Metabolite",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4], [scatter5, line5], [scatter6, line6]],
    ["g6p", "f6p", "fdp", "dhap", "g3p", "dpg13"],
)
f

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Metabolites_NADH_NAD_abs.pdf"), f)

#=
Plot lower glycolysis metabolites
=#
x = exp.(nadh)./exp.(nad)
a = exp.([pg3]...)
b = exp.([pg2]...)
c = exp.([pep]...)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "Metabolite concentration [M]"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="Metabolite",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3]],
    ["p3g", "pg2", "pep"],
)
f

FileIO.save(joinpath(imgpath, "Lower_Glycolysis_Metabolites_NADH_NAD_abs.pdf"), f)

#=
Plot upper glycolysis proteome
=#
x = exp.(nadh)./exp.(nad)
a = (pts)
b = (pgi)
c = (pfk)
d = (fba)
e = (tpi)
g = (gapd)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d)
line4 = lines!(ax, x, d)
scatter5 = scatter!(ax, x, e)
line5 = lines!(ax, x, e)
scatter6 = scatter!(ax, x, g)
line6 = lines!(ax, x, g)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "Enzyme concentration [g enz / g DW]"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="Enzyme",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4], [scatter5, line5], [scatter6, line6]],
    ["pts", "pgi", "pfk", "fba", "tpi", "gapd"],
)
f

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Proteome_NADH_NAD_abs.pdf"), f)

#=
Plot lower glycolysis proteome
=#
x = exp.(nadh)./exp.(nad)
a = (pgk)
b = (pgm)
c = (pyk)
d = (eno)
e = (ldh)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d)
line4 = lines!(ax, x, d)
scatter5 = scatter!(ax, x, e)
line5 = lines!(ax, x, e)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "Enzyme concentration [g enz / g DW]"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="Enzyme",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4], [scatter5, line5]],
    ["pgk", "pgm", "pyk", "eno", "ldh"],
)
f

FileIO.save(joinpath(imgpath, "Lower_Glycolysis_Proteome_NADH_NAD_abs.pdf"), f)

#=
Plot Fluxes of reactions in upper glycolysis
=#
x = exp.(nadh)./exp.(nad)
a = (v_pts)
b = (v_pgi)
c = (v_pfk)
d = (v_fba)
e = (v_tpi)
g = (v_gapd)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d)
line4 = lines!(ax, x, d)
scatter5 = scatter!(ax, x, e)
line5 = lines!(ax, x, e)
scatter6 = scatter!(ax, x, g)
line6 = lines!(ax, x, g)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "Fluxes [mmol/gDW/h]"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="Flux",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4], [scatter5, line5], [scatter6, line6]],
    ["v_pts", "v_pgi", "v_pfk", "v_fba", "v_tpi", "v_gapd"],
)
f

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Fluxes_NADH_NAD_abs.pdf"), f)

#=
Plot dG of reactions in upper glycolysis
=#
x = exp.(nadh)./exp.(nad)
a = (dg_pts)
b = (dg_pgi)
c = (dg_pfk)
d = (dg_fba)
e = (dg_tpi)
g = (dg_gapd)

f = Figure()
ax = Axis(f[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c)
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d)
line4 = lines!(ax, x, d)
scatter15 = scatter!(ax, x, e)
line5 = lines!(ax, x, e)
scatter6 = scatter!(ax, x, g)
line6 = lines!(ax, x, g)

ax.xlabel = "NADH/NAD ratio"
ax.ylabel = "dG of reaction"
xscale = log10
f[1, 2] = Legend(
    f,
    legendlabel="dG of reaction",
    [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4],[scatter5, line5], [scatter6, line6]],
    ["dg_pts", "dg_pgi", "dg_pfk", "dg_fba", "dg_tpi", "dg_gapd"],
)
f

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_dG_NADH_NAD_abs.pdf"), f)