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
ax.ylabel = "Proteome fraction"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glycolysis_Glc.pdf"), f)

## Glutamine synthesis enzymes ##

fracs = [a2, b2, c2, d2, e2, f2, g2, h2, i2]
fraclabels = ["pgi", "pfk", "fba", "tpi", "gapd", "pgk", "pgm", "eno", "pyk"]

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
ax.ylabel = "Proteome fraction"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Glutamine_syn_Glc.pdf"), f)

## Import + Maintanance enzymes ##

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
ax.ylabel = "Proteome fraction"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Proteome_Import_Mainanance_Glc.pdf"), f)
