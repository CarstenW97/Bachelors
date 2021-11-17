using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("docs", "imgs", "EMP-model")
resultspath = joinpath("docs", "results", "EMP-model")

## PTS concentration - Glucose concentration ##

x = [exp(results_emp["glc_ext"]) for results_emp in res_emp]
y = [exp(results_emp["pts"]) for results_emp in res_emp]

fig = Figure()
ax = Axis(fig[1, 1])
scatter_glc = scatter!(ax, x, y)
line_glc = lines!(ax, x, y)

ax.xlabel = "External glucose concetration [mM]"
ax.ylabel = "PTS concentration [g enz / g DW]"

fig

FileIO.save(joinpath(imgpath, "EMP_PTScon_Glc.pdf"), fig)

## Lactate flux - Glucose concentration ##

x = [exp(results_emp["glc_ext"]) for results_emp in res_emp]
y = [results_emp["v_ldh"] for results_emp in res_emp]

fig = Figure()
ax = Axis(fig[1, 1])
scatter_glc = scatter!(ax, x, y)
line_glc = lines!(ax, x, y)

ax.xlabel = "External glucose concetration [mM]"
ax.ylabel = "Lactate production flux [mmol/gDW/h]"

fig

FileIO.save(joinpath(imgpath, "EMP_Lacflux_Glc.pdf"), fig)

## Glutamine flux - Glucose concentration ##

x = [exp(results_emp["glc_ext"]) for results_emp in res_emp]
y = [results_emp["v_glnsyn"] for results_emp in res_emp]

fig = Figure()
ax = Axis(fig[1, 1])
scatter_glc = scatter!(ax, x, y)
line_glc = lines!(ax, x, y)

ax.xlabel = "External glucose concetration [mM]"
ax.ylabel = "Glutamine production flux [mmol/gDW/h]"

fig

FileIO.save(joinpath(imgpath, "EMP_Glnflux_Glc.pdf"), fig)

## Proteome - glucose concentration ##

x = [exp(results_ed["glc_ext"]) for results_ed in res_emp]
a = [results_ed["pts"] for results_ed in res_emp]
b = [results_ed["emp"] for results_ed in res_emp]
c = [results_ed["pyk"] for results_ed in res_emp]
d = [results_ed["ldh"] for results_ed in res_emp]
e = [results_ed["ppc"] for results_ed in res_emp]
f = [results_ed["akgsyn"] for results_ed in res_emp]
g = [results_ed["gdhm"] for results_ed in res_emp]
h = [results_ed["burn"] for results_ed in res_emp]
i = [results_ed["nadtrdh"] for results_ed in res_emp]
j = [results_ed["lp"] for results_ed in res_emp]
k = [results_ed["glnsyn"] for results_ed in res_emp]
l = [results_ed["nh3_diff"] for results_ed in res_emp]
m = [results_ed["co2_diff"] for results_ed in res_emp]

fracs = [a, b, c, d, e, f, g, h, i, j, k, l, m]
fraclabels = [
    "PTS",
    "EMP",
    "Pyruvate kinase",
    "Lactate dehydrogenase",
    "Phosphoenolpyruvate carboxylase",
    "alpah-Ketogluterate synthesis enzyme",
    "Glutamate synthesis enzyme",
    "ATP burning enzyme",
    "Lactate permease",
    "Glutamine synthesis enzyme",
    "NH3 permease",
    "CO2 permease",
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

FileIO.save(joinpath(imgpath, "EMP_Proteome_Glc.pdf"), f)
