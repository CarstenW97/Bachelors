using JuMP
using JSON

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel

using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("docs", "imgs", "EMP-model")
resultspath = joinpath("docs", "results", "EMP-model")

## Proteome - glucose concentration ##

x  = [] # glc_ext

a1 = [] # pts 
b1 = [] # pgi 
c1 = [] # pfk 
d1 = [] # fba 
e1 = [] # tpi 
f1 = [] # gapd
g1 = [] # pgk 
h1 = [] # pgm 
i1 = [] # eno 
j1 = [] # pyk 

a2 = [] # ppc
b2 = [] # pdh
c2 = [] # pfl
d2 = [] # cs
e2 = [] # aconta
f2 = [] # acontb
g2 = [] # icdh
h2 = [] # gludy
i2 = [] # glns

c3 = [] # ldh
d3 = [] # lact
e3 = [] # fort
f3 = [] # ptar
g3 = [] # ackr
h3 = [] # act
i3 = [] # acald
j3 = [] # alcd
k3 = [] # etoht
l3 = [] # nh4t

for var in range(0.005, 0.5; length=10)

    m = GlnModel.gln_model(
        glc_ext = var,
        lac_ext = 1e-2,
        nh4_ext = 0.01,
        ac_ext = 1e-16,
        etoh_ext = 1e-6,
        atp_adp_ratio = 10.0,
        nadh_nad_ratio = 0.2,
    )
    
    push!(x, value(var))

    push!(a1, value(m[:pts]))
    push!(b1, value(m[:pgi]))
    push!(c1, value(m[:pfk]))
    push!(d1, value(m[:fba]))
    push!(e1, value(m[:tpi]))
    push!(f1, value(m[:gapd]))
    push!(g1, value(m[:pgk]))
    push!(h1, value(m[:pgm]))
    push!(i1, value(m[:eno]))
    push!(j1, value(m[:pyk]))

    push!(a2, value(m[:ppc]))
    push!(b2, value(m[:pdh]))
    push!(c2, value(m[:pfl]))
    push!(d2, value(m[:cs]))
    push!(e2, value(m[:aconta]))
    push!(f2, value(m[:acontb]))
    push!(g2, value(m[:icdh]))
    push!(h2, value(m[:gludy]))
    push!(i2, value(m[:glns]))

    push!(c3, value(m[:ldh]))
    push!(d3, value(m[:lact]))
    push!(e3, value(m[:fort]))
    push!(f3, value(m[:ptar]))
    push!(g3, value(m[:ackr]))
    push!(h3, value(m[:act]))
    push!(i3, value(m[:acald]))
    push!(j3, value(m[:alcd]))
    push!(k3, value(m[:etoht]))
    push!(l3, value(m[:nh4t]))

end

## Glycolysis enzymes ##

fracs = [a1, b1, c1, d1, e1, f1, g1, h1, i1, j1]
fraclabels = ["pts", "pgi", "pfk", "fba", "tpi", "gapd", "pgk", "pgm", "eno", "pyk"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1]);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "Glucose concentration"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glycolysis_Glc.pdf"), f)

## Glutamine synthesis enzymes ##

fracs = [a2, b2, c2, d2, e2, f2, g2, h2, i2]
fraclabels = ["ppc", "pdh", "pfl", "cs", "aconta", "acontb", "icdh", "gludy", "glns"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1]);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "Glucose concentration"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glutamine_syn_Glc.pdf"), f)

## Transport + Maintanance enzymes ##

fracs = [c3, d3, e3, f3, g3, h3, i3, j3, k3, l3]
fraclabels = ["ldh", "lact", "fort", "ptar", "ackr", "act", "acald", "alcd", "etoht", "nh4t"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1]);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "Glucose concentration"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Transport_Mainanance_Glc.pdf"), f)





## ATP/ADP-Ratio ##

x  = [] # atp_adp_ratio

    # Enzymes
a1 = [] # pts 
b1 = [] # pgi 
c1 = [] # pfk 
d1 = [] # fba 
e1 = [] # tpi 
f1 = [] # gapd
g1 = [] # pgk 
h1 = [] # pgm 
i1 = [] # eno 
j1 = [] # pyk 

a2 = [] # ppc
b2 = [] # pdh
c2 = [] # pfl
d2 = [] # cs
e2 = [] # aconta
f2 = [] # acontb
g2 = [] # icdh
h2 = [] # gludy
i2 = [] # glns

c3 = [] # ldh
d3 = [] # lact
e3 = [] # fort
f3 = [] # ptar
g3 = [] # ackr
h3 = [] # act
i3 = [] # acald
j3 = [] # alcd
k3 = [] # etoht
l3 = [] # nh4t

    # Metabolites
a4 = [] # g6p
b4 = [] # f6p
c4 = [] # fdp
d4 = [] # dhap
e4 = [] # g3p
f4 = [] # dpg13
g4 = [] # pg3
h4 = [] # pg2
i4 = [] # pep
j4 = [] # pyr

a5 = [] # lac
b5 = [] # accoa
c5 = [] # oaa
d5 = [] # cit
e5 = [] # acon
f5 = [] # icit
g5 = [] # akg
h5 = [] # nh4
i5 = [] # glu
j5 = [] # gln

a6 = [] # adp 
b6 = [] # atp 
c6 = [] # nad 
d6 = [] # nadh
e6 = [] # acetald
f6 = [] # etoh
g6 = [] # actp
h6 = [] # ac
i6 = [] # formate

for var in range(0.01, 10.0; length=100)

    m = GlnModel.gln_model(
        glc_ext = 0.05,
        lac_ext = 1e-2,
        nh4_ext = 0.01,
        ac_ext = 1e-16,
        etoh_ext = 1e-6,
        atp_adp_ratio = var,
        nadh_nad_ratio = 0.2,
    )
    
    push!(x, value(var))

        # Enzymes
    push!(a1, value(m[:pts]))
    push!(b1, value(m[:pgi]))
    push!(c1, value(m[:pfk]))
    push!(d1, value(m[:fba]))
    push!(e1, value(m[:tpi]))
    push!(f1, value(m[:gapd]))
    push!(g1, value(m[:pgk]))
    push!(h1, value(m[:pgm]))
    push!(i1, value(m[:eno]))
    push!(j1, value(m[:pyk]))

    push!(a2, value(m[:ppc]))
    push!(b2, value(m[:pdh]))
    push!(c2, value(m[:pfl]))
    push!(d2, value(m[:cs]))
    push!(e2, value(m[:aconta]))
    push!(f2, value(m[:acontb]))
    push!(g2, value(m[:icdh]))
    push!(h2, value(m[:gludy]))
    push!(i2, value(m[:glns]))

    push!(c3, value(m[:ldh]))
    push!(d3, value(m[:lact]))
    push!(e3, value(m[:fort]))
    push!(f3, value(m[:ptar]))
    push!(g3, value(m[:ackr]))
    push!(h3, value(m[:act]))
    push!(i3, value(m[:acald]))
    push!(j3, value(m[:alcd]))
    push!(k3, value(m[:etoht]))
    push!(l3, value(m[:nh4t]))

        # Metaolites
    push!(a4, value(m[:g6p]))
    push!(b4, value(m[:f6p]))
    push!(c4, value(m[:fdp]))
    push!(d4, value(m[:dhap]))
    push!(e4, value(m[:g3p]))
    push!(f4, value(m[:dpg13]))
    push!(g4, value(m[:pg3]))
    push!(h4, value(m[:pg2]))
    push!(i4, value(m[:pep]))
    push!(j4, value(m[:pyr]))
    
    push!(a5, value(m[:lac]))
    push!(b5, value(m[:accoa]))
    push!(c5, value(m[:oaa]))
    push!(d5, value(m[:cit]))
    push!(e5, value(m[:acon]))
    push!(f5, value(m[:icit]))
    push!(g5, value(m[:akg]))
    push!(h5, value(m[:nh4]))
    push!(i5, value(m[:glu]))
    push!(j5, value(m[:gln]))
    
    push!(a6, value(m[:adp])) 
    push!(b6, value(m[:atp])) 
    push!(c6, value(m[:nad])) 
    push!(d6, value(m[:nadh]))
    push!(e6, value(m[:acetald]))
    push!(f6, value(m[:etoh]))
    push!(g6, value(m[:actp]))
    push!(h6, value(m[:ac]))
    push!(i6, value(m[:formate]))

end

## Glycolysis enzymes ##

fracs = [a1, b1, c1, d1, e1, f1, g1, h1, i1, j1]
fraclabels = ["pts", "pgi", "pfk", "fba", "tpi", "gapd", "pgk", "pgm", "eno", "pyk"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glycolysis_ATP_ADP_Ratio_log.pdf"), f)

## Glutamine synthesis enzymes ##

fracs = [a2, b2, c2, d2, e2, f2, g2, h2, i2]
fraclabels = ["ppc", "pdh", "pfl", "cs", "aconta", "acontb", "icdh", "gludy", "glns"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glutamine_syn_ATP_ADP_Ratio_log.pdf"), f)

## Transport + Maintanance enzymes ##

fracs = [c3, d3, e3, f3, g3, h3, i3, j3, k3, l3]
fraclabels = ["ldh", "lact", "fort", "ptar", "ackr", "act", "acald", "alcd", "etoht", "nh4t"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Transport_Mainanance_ATP_ADP_Ratio_log.pdf"), f)

# Metabolite plots #

fracs = [a4, b4, c4, d4, e4, f4, g4, h4, i4, j4]
fraclabels = ["g6p", "f6p", "fdp", "dhap", "g3p", "dpg13", "pg3", "pg2", "pep", "pyr"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Metabolite concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Metabolite_1_ATP_ADP_Ratio_log.pdf"), f)

#########

fracs = [a5, b5, c5, d5, e5, f5, g5, h5, i5, j5]
fraclabels = ["lac", "accoa", "oaa", "cit", "acon", "icit", "akg", "nh4", "glu", "gln"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Metabolite concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Metabolite_2_ATP_ADP_Ratio_log.pdf"), f)

#########

fracs = [a6, b6, c6, d6, e6, f6, g6, h6, i6]
fraclabels = ["adp", "atp", "nad", "nadh", "acetald", "etoh", "actp", "ac", "formate"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Metabolite concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Metabolite_3_ATP_ADP_Ratio_log.pdf"), f)





## NAD/NADH-Ratio ##

x  = [] # nadh_nad_ratio

    # Enzymes
a1 = [] # pts 
b1 = [] # pgi 
c1 = [] # pfk 
d1 = [] # fba 
e1 = [] # tpi 
f1 = [] # gapd
g1 = [] # pgk 
h1 = [] # pgm 
i1 = [] # eno 
j1 = [] # pyk 

a2 = [] # ppc
b2 = [] # pdh
c2 = [] # pfl
d2 = [] # cs
e2 = [] # aconta
f2 = [] # acontb
g2 = [] # icdh
h2 = [] # gludy
i2 = [] # glns

c3 = [] # ldh
d3 = [] # lact
e3 = [] # fort
f3 = [] # ptar
g3 = [] # ackr
h3 = [] # act
i3 = [] # acald
j3 = [] # alcd
k3 = [] # etoht
l3 = [] # nh4t

    # Metabolites
a4 = [] # g6p
b4 = [] # f6p
c4 = [] # fdp
d4 = [] # dhap
e4 = [] # g3p
f4 = [] # dpg13
g4 = [] # pg3
h4 = [] # pg2
i4 = [] # pep
j4 = [] # pyr

a5 = [] # lac
b5 = [] # accoa
c5 = [] # oaa
d5 = [] # cit
e5 = [] # acon
f5 = [] # icit
g5 = [] # akg
h5 = [] # nh4
i5 = [] # glu
j5 = [] # gln

a6 = [] # adp 
b6 = [] # atp 
c6 = [] # nad 
d6 = [] # nadh
e6 = [] # acetald
f6 = [] # etoh
g6 = [] # actp
h6 = [] # ac
i6 = [] # formate

for var in range(0.01, 10.0; length=100)

    m = GlnModel.gln_model(
        glc_ext = 0.05,
        lac_ext = 1e-2,
        nh4_ext = 0.01,
        ac_ext = 1e-16,
        etoh_ext = 1e-6,
        atp_adp_ratio = 10.0,
        nadh_nad_ratio = var,
    )
    
    push!(x, value(var))

        # Enzymes
    push!(a1, value(m[:pts]))
    push!(b1, value(m[:pgi]))
    push!(c1, value(m[:pfk]))
    push!(d1, value(m[:fba]))
    push!(e1, value(m[:tpi]))
    push!(f1, value(m[:gapd]))
    push!(g1, value(m[:pgk]))
    push!(h1, value(m[:pgm]))
    push!(i1, value(m[:eno]))
    push!(j1, value(m[:pyk]))

    push!(a2, value(m[:ppc]))
    push!(b2, value(m[:pdh]))
    push!(c2, value(m[:pfl]))
    push!(d2, value(m[:cs]))
    push!(e2, value(m[:aconta]))
    push!(f2, value(m[:acontb]))
    push!(g2, value(m[:icdh]))
    push!(h2, value(m[:gludy]))
    push!(i2, value(m[:glns]))

    push!(c3, value(m[:ldh]))
    push!(d3, value(m[:lact]))
    push!(e3, value(m[:fort]))
    push!(f3, value(m[:ptar]))
    push!(g3, value(m[:ackr]))
    push!(h3, value(m[:act]))
    push!(i3, value(m[:acald]))
    push!(j3, value(m[:alcd]))
    push!(k3, value(m[:etoht]))
    push!(l3, value(m[:nh4t]))

        # Metaolites
    push!(a4, value(m[:g6p]))
    push!(b4, value(m[:f6p]))
    push!(c4, value(m[:fdp]))
    push!(d4, value(m[:dhap]))
    push!(e4, value(m[:g3p]))
    push!(f4, value(m[:dpg13]))
    push!(g4, value(m[:pg3]))
    push!(h4, value(m[:pg2]))
    push!(i4, value(m[:pep]))
    push!(j4, value(m[:pyr]))
    
    push!(a5, value(m[:lac]))
    push!(b5, value(m[:accoa]))
    push!(c5, value(m[:oaa]))
    push!(d5, value(m[:cit]))
    push!(e5, value(m[:acon]))
    push!(f5, value(m[:icit]))
    push!(g5, value(m[:akg]))
    push!(h5, value(m[:nh4]))
    push!(i5, value(m[:glu]))
    push!(j5, value(m[:gln]))
    
    push!(a6, value(m[:adp])) 
    push!(b6, value(m[:atp])) 
    push!(c6, value(m[:nad])) 
    push!(d6, value(m[:nadh]))
    push!(e6, value(m[:acetald]))
    push!(f6, value(m[:etoh]))
    push!(g6, value(m[:actp]))
    push!(h6, value(m[:ac]))
    push!(i6, value(m[:formate]))

end

## Glycolysis enzymes ##

fracs = [a1, b1, c1, d1, e1, f1, g1, h1, i1, j1]
fraclabels = ["pts", "pgi", "pfk", "fba", "tpi", "gapd", "pgk", "pgm", "eno", "pyk"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "NAD/NADH-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glycolysis_NAD_NADH_Ratio_log.pdf"), f)

## Glutamine synthesis enzymes ##

fracs = [a2, b2, c2, d2, e2, f2, g2, h2, i2]
fraclabels = ["ppc", "pdh", "pfl", "cs", "aconta", "acontb", "icdh", "gludy", "glns"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "NAD/NADH-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glutamine_syn_NAD_NADH_Ratio_log.pdf"), f)

## Transport + Maintanance enzymes  #full# ##

fracs = [c3, d3, e3, f3, g3, h3, i3, j3, k3, l3]
fraclabels = ["ldh", "lact", "fort", "ptar", "ackr", "act", "acald", "alcd", "etoht", "nh4t"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "NAD/NADH-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Transport_Mainanance_NAD_NADH_Ratio_log.pdf"), f)

## Transport + Maintanance enzymes  #1# ##

fracs = [c3, d3, e3, f3, g3]
fraclabels = ["ldh", "lact", "fort", "ptar", "ackr"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "NAD/NADH-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Transport_Mainanance_NAD_NADH_Ratio_log_1.pdf"), f)

## Transport + Maintanance enzymes  #2# ##

fracs = [h3, i3, j3, k3, l3]
fraclabels = ["act", "acald", "alcd", "etoht", "nh4t"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "NAD/NADH-Ratio"
ax.ylabel = "Relative Percental Enzyme concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Transport_Mainanance_NAD_NADH_Ratio_log_2.pdf"), f)


# Metabolite plots #

fracs = [a4, b4, c4, d4, e4, f4, g4, h4, i4, j4]
fraclabels = ["g6p", "f6p", "fdp", "dhap", "g3p", "dpg13", "pg3", "pg2", "pep", "pyr"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Metabolite concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Metabolite_1_NAD_NADH_Ratio_log.pdf"), f)

#########

fracs = [a5, b5, c5, d5, e5, f5, g5, h5, i5, j5]
fraclabels = ["lac", "accoa", "oaa", "cit", "acon", "icit", "akg", "nh4", "glu", "gln"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Metabolite concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Metabolite_2_NAD_NADH_Ratio_log.pdf"), f)

#########

fracs = [a6, b6, c6, d6, e6, f6, g6, h6, i6]
fraclabels = ["adp", "atp", "nad", "nadh", "acetald", "etoh", "actp", "ac", "formate"]

function plot_proteome(ax, x, fracs, fraclabels)
    z = cumsum(hcat(fracs...)'; dims = 1)
    z = [zeros(size(z, 2))'; z ./ z[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, z[i, :], z[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1], xscale=log10);
plot_proteome(ax, x, fracs, fraclabels)
ax.xlabel = "ATP/ADP-Ratio"
ax.ylabel = "Relative Percental Metabolite concentration"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Metabolite_3_NAD_NADH_Ratio_log.pdf"), f)