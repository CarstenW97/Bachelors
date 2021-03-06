using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("docs", "imgs", "ED-model")
resultspath = joinpath("docs", "results", "ED-model")

x = [results_ed["glc_ext"] for results_ed in res_ed]
y = [results_ed["mu"] for results_ed in res_ed]

fig = Figure()
ax = Axis(fig[1, 1])
scatter_glc = scatter!(ax, x, y)
line_glc = lines!(ax, x, y)

ax.xlabel = "External glucose concetration [mM]"
ax.ylabel = "Biomass [1/h]"

fig

FileIO.save(joinpath(imgpath, "Glu_ED_Growth.pdf"), fig)

## ##

x = [results_ed["mu"] for results_ed in res_ed]
a = [results_ed["g6p"] for results_ed in res_ed]
b = [results_ed["glu"] for results_ed in res_ed]

fig = Figure()
ax = Axis(fig[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)

ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Metabolite concentration"

fig[1, 1] = Legend(
    fig,
    [[scatter1, line1], [scatter2, line2]],
    ["Glucose-6-Phosphate", "Glutamate"],
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :bottom,
    labelsize = 14,
)
fig

FileIO.save(joinpath(imgpath, "Glu_ED_G6P+GLU.pdf"), fig)

## Enzyme concentration band plot ##

x = [results_ed["mu"] for results_ed in res_ed]
c = [results_ed["pts"] for results_ed in res_ed]
d = [results_ed["ed"] for results_ed in res_ed]
e = [results_ed["pyk"] for results_ed in res_ed]
f = [results_ed["ldh"] for results_ed in res_ed]
g = [results_ed["burn"] for results_ed in res_ed]
h = [results_ed["lp"] for results_ed in res_ed]

fracs = [c, d, e, f, g, h]
fraclabels = [
    "PTS",
    "ED",
    "Pyruvate kinase",
    "Lactate dehydrogenase",
    "ATP burning enzyme",
    "Lactate permease",
]
function plot_proteome(ax, x, fracs, fraclabels)
    d = cumsum(hcat(fracs...)'; dims = 1)
    d = [zeros(size(d, 2))'; d ./ d[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, d[i, :], d[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1]);
plot_proteome(ax, 1:10, fracs, fraclabels)
ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Proteome fraction"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Glu_ED_Enzymeband.pdf"), f)

## Rates ##

x = [results_ed["mu"] for results_ed in res_ed]
i = [results_ed["v_pts"] for results_ed in res_ed]
j = [results_ed["v_pyk"] for results_ed in res_ed]
k = [results_ed["v_ldh"] for results_ed in res_ed]
l = [results_ed["v_burn"] for results_ed in res_ed]
m = [results_ed["v_lp"] for results_ed in res_ed]

fig = Figure()
ax = Axis(fig[1, 1])
scatter1 = scatter!(ax, x, i)
line1 = lines!(ax, x, i)
scatter2 = scatter!(ax, x, j)
line2 = lines!(ax, x, j)
scatter3 = scatter!(ax, x, k)
line3 = lines!(ax, x, k)
scatter4 = scatter!(ax, x, l)
line4 = lines!(ax, x, l)
scatter5 = scatter!(ax, x, m)
line5 = lines!(ax, x, m)

ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Fluxes [mmol/gDW/h]"

fig[1, 1] = Legend(
    fig,
    [
        [scatter1, line1],
        [scatter2, line2],
        [scatter3, line3],
        [scatter4, line4],
        [scatter5, line5],
    ],
    ["v_pts", "v_pyk", "v_ldh", "v_burn", "v_lp"],
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :bottom,
    labelsize = 14,
)
fig

FileIO.save(joinpath(imgpath, "Glu_ED_Rates.pdf"), fig)
