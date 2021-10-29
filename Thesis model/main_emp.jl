include(joinpath("Thesis model", "glu_emp_model.jl"))
import .Glu_EMP_Model

res_emp = []
for var in [0.0001:0.0005:0.005;]
  push!(res_emp, Glu_EMP_Model.glu_emp_model(var))
end

################################################################

using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Thesis model", "Figures" , "EMP_Model")
resultspath = joinpath("Thesis model", "Results", "EMP_Model")

res_df = Array{String, Float64}(res_emp)
df = DataFrame[]
for i in 1:100
  push!(df, res_emp[i])
end
CSV.write(joinpath(resultspath, "Results_EMP.csv"), res_emp[1])


x = [results_emp["mu"] for results_emp in res_emp]
a = [results_emp["g6p"] for results_emp in res_emp]
b = [results_emp["glu"] for results_emp in res_emp]

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

FileIO.save(joinpath(imgpath, "Glu_EMP.pdf"), fig) 

## Enzyme concentration band plot ##

x = [results_emp["mu"] for results_emp in res_emp]
c = [results_emp["pts"] for results_emp in res_emp]
d = [results_emp["emp"] for results_emp in res_emp]
e = [results_emp["pyk"] for results_emp in res_emp]
f = [results_emp["ldh"] for results_emp in res_emp]
g = [results_emp["burn"] for results_emp in res_emp]
h = [results_emp["lp"] for results_emp in res_emp]

fracs = [c, d, e, f, g, h]
fraclabels = ["PTS", "EMP", "Pyruvate kinase", "Lactate dehydrogenase", 
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

FileIO.save(joinpath(imgpath, "Enzyme_band.pdf"), f)


