using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("docs", "imgs", "EMP-model")
resultspath = joinpath("docs", "results", "EMP-model")

## Proteome - glucose concentration ##

x = value(m[glc_e])
a1 = value(m[:pts])
b1 = value(m[:pgi])
c1 = value(m[:pfk])
d1 = value(m[:fba])
e1 = value(m[:tpi])
f1 = value(m[:gapd])
g1 = value(m[:pgk])
h1 = value(m[:pgm])
i1 = value(m[:eno])
j1 = value(m[:pyk])
k1 = value(m[:ppc])
l1 = value(m[:ldh])
m1 = value(m[:pdh])
n1 = value(m[:lact])
o1 = value(m[:cs])
p1 = value(m[:aconta])
q1 = value(m[:acontb])
r1 = value(m[:icdh])
s1 = value(m[:gludy])
t1 = value(m[:glns])
u1 = value(m[:nh4t])
v1 = value(m[:acald])
w1 = value(m[:alcd])
a2 = value(m[:ethot])
b2 = value(m[:ackr])
c2 = value(m[:ptar])
d2 = value(m[:act])
e2 = value(m[:pfl])
f2 = value(m[:fort])

fracs = [
    a1,
    b1,
    c1,
    d1,
    e1,
    f1,
    g1,
    h1,
    i1,
    j1,
    k1,
    l1,
    m1,
    n1,
    o1,
    p1,
    q1,
    r1,
    s1,
    t1,
    u1,
    v1,
    w1,
    a2,
    b2,
    c2,
    d2,
    e2,
    f2,
]
fraclabels = [
    "pts",
    "pgi",
    "pfk",
    "fba",
    "tpi",
    "gapd",
    "pgk",
    "pgm",
    "eno",
    "pyk",
    "ppc",
    "ldh",
    "pdh",
    "lact",
    "cs",
    "aconta",
    "acontb",
    "icdh",
    "gludy",
    "glns",
    "nh4t",
    "acald",
    "alcd",
    "etoht",
    "ackr",
    "ptar",
    "act",
    "pfl",
    "fort",
]
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

FileIO.save(joinpath(imgpath, "Proteome_Glc.pdf"), f)
