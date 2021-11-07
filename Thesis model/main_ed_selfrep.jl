include(joinpath("Thesis model", "glu_ed_selfrep_model.jl"))
import .Glu_ED_Selfrep_Model

res_ed_selfrep = []
for var in [0.0001:0.0005:0.005;]
  push!(res_ed_selfrep, Glu_ED_Selfrep_Model.glu_ed_selfrep_model(var))
end

################################################################

using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Thesis model", "Figures" , "ED_Model")
resultspath = joinpath("Thesis model", "Results", "ED_Model")

#CSV.write(joinpath(resultspath, "Results_ED_Selfrep.csv"), res_ed[1])

x = [results_ed_selfrep["glc_ext"] for results_ed_selfrep in res_ed_selfrep]
y = [results_ed_selfrep["mu"] for results_ed_selfrep in res_ed_selfrep]

fig = Figure()
ax = Axis(fig[1,1])
scatter_glc = scatter!(ax, x, y)
line_glc = lines!(ax, x, y)

ax.xlabel = "External glucose concetration [mM]"
ax.ylabel = "Biomass [1/h]"

fig

FileIO.save(joinpath(imgpath, "Glu_ED_Selfrep_Growth.pdf"), fig)

## ##

x = [results_ed_selfrep["mu"]  for results_ed_selfrep in res_ed_selfrep]
a = [results_ed_selfrep["g6p"] for results_ed_selfrep in res_ed_selfrep]
b = [results_ed_selfrep["glu"] for results_ed_selfrep in res_ed_selfrep]

fig = Figure()
ax = Axis(fig[1,1])
scatter1 = scatter!(ax, x, a)
line1 =lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 =lines!(ax, x, b)

ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Metabolite concentration"

fig[1, 1] = Legend(fig, [[scatter1, line1], [scatter2, line2]], ["Glucose-6-Phosphate", "Glutamate"],
    tellheight = false, tellwidth = false, halign = :left, valign = :bottom, labelsize = 14)
fig

FileIO.save(joinpath(imgpath, "Glu_ED_G6P+GLU_Selfrep.pdf"), fig) 

## Enzyme concentration band plot ##

x = [results_ed_selfrep["mu"]   for results_ed_selfrep in res_ed_selfrep]
c = [results_ed_selfrep["pts"]  for results_ed_selfrep in res_ed_selfrep]
d = [results_ed_selfrep["ed"]  for results_ed_selfrep in res_ed_selfrep]
e = [results_ed_selfrep["pyk"]  for results_ed_selfrep in res_ed_selfrep]
f = [results_ed_selfrep["ldh"]  for results_ed_selfrep in res_ed_selfrep]
g = [results_ed_selfrep["burn"] for results_ed_selfrep in res_ed_selfrep]
h = [results_ed_selfrep["lp"]   for results_ed_selfrep in res_ed_selfrep]

fracs = [c, d, e, f, g, h]
fraclabels = ["PTS", "ED", "Pyruvate kinase", "Lactate dehydrogenase", 
              "ATP burning enzyme", "Lactate permease"]
function plot_proteome(ax, x, fracs, fraclabels)
    d = cumsum(hcat(fracs...)'; dims=1)
    d = [zeros(size(d, 2))'; d ./ d[end,:]']
    for i in 1:length(fraclabels)
        band!(ax, x, d[i,:], d[i+1,:], label=fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1,1]);
plot_proteome(ax, 1:10, fracs, fraclabels)
ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Proteome fraction"
f[1,2] = Legend(f, ax, "Protein",
                    unique=true,  framevisible = false)
f

FileIO.save(joinpath(imgpath, "Glu_ED_Selfrep_Enzymeband.pdf"), f)

## Rates ##

x = [results_ed_selfrep["mu"]     for results_em_selfrep in res_ed_selfrep]
i = [results_ed_selfrep["v_pts"]  for results_em_selfrep in res_ed_selfrep]
j = [results_ed_selfrep["v_pyk"]  for results_em_selfrep in res_ed_selfrep]
k = [results_ed_selfrep["v_ldh"]  for results_em_selfrep in res_ed_selfrep]
l = [results_ed_selfrep["v_burn"] for results_em_selfrep in res_ed_selfrep]
m = [results_ed_selfrep["v_lp"]   for results_em_selfrep in res_ed_selfrep]

fig = Figure()
ax = Axis(fig[1,1])
scatter1 = scatter!(ax, x, i)
line1 =lines!(ax, x, i)
scatter2 = scatter!(ax, x, j)
line2 =lines!(ax, x, j)
scatter3 = scatter!(ax, x, k)
line3 =lines!(ax, x, k)
scatter4 = scatter!(ax, x, l)
line4 =lines!(ax, x, l)
scatter5 = scatter!(ax, x, m)
line5 =lines!(ax, x, m)

ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Fluxes [mmol/gDW/h]"

fig[1, 1] = Legend(fig, [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4], [scatter5, line5]], 
                         ["v_pts", "v_pyk", "v_ldh", "v_burn", "v_lp"],
    tellheight = false, tellwidth = false, halign = :left, valign = :bottom, labelsize = 14)
fig

FileIO.save(joinpath(imgpath, "Glu_ED_Selfrep_Rates.pdf"), fig)