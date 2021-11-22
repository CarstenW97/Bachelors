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
@showprogress for var in range(log(ratio_lb); stop=log(ratio_ub), length=20)
    nadh_nad_ratio = exp(var)

    gln_model = GlnModel.gln_model(
        glc_ext = 50e-3,
        lac_ext = 1e-4,
        nh4_ext = 10e-3,
        ac_ext = 1e-16,
        etoh_ext = 1e-4,
        # atp_adp_ratio = -1, # unconstrained
        nadh_nad_ratio = nadh_nad_ratio,
        num_ms = 100,
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
    ylabel="Relative Gibbs dissipation",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Upper_Glycolysis_Fluxes_NADH_NAD.pdf"), f)

###############

fracs = hcat([pts, pgi, pfk, fba, tpi, gapd]...)

fracs = hcat([pts, pgi, pfk, fba, tpi, gapd]...)
fraclabels = ["pts", "pgi", "pfk", "fba", "tpi", "gapd"]
x = exp.(nadh)./exp.(nad)

f = plot_proteome(x, fracs, fraclabels;
    legendlabel="Enzyme",
    xlabel="NADH/NAD ratio",
    ylabel="Relative Enzyme concentration",
    xscale=log10,
)

FileIO.save(joinpath(imgpath, "Test_NADH_NAD.pdf"), f)