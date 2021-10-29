include(joinpath("Thesis model", "glu_emp_model.jl"))
import .Glu_ED_Model

res_ed = []
for var in [20]
  push!(res_ed, Glu_ED_Model.glu_ed_model(var))
end

################################################################

using CSV
CSV.write(joinpath(resultpath, "results.csv") results_ed)

################################################################

using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Thesis model", "Figures" , "ED_Model")
resultspath = joinpath("Thesis model", "Results", "ED_Model")

x = [results_ed["mu"] for results_ed in res_ed]
a = [results_ed["g6p"] for results_ed in res_ed]
b = [results_ed["glu"] for results_ed in res_ed]

a = Float64.(a)
b = Float64.(b)

fig = Figure()
ax = Axis(fig[1,1])
scatter1 = scatter!(ax, x, a)
line1 =lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 =lines!(ax, x, b)

ax.xlabel = "Biomass [1/h]"
ax.ylabel = "Metabolite concentration"

fig[1, 1] = Legend(fig, [[scatter1, line1], [scatter2, line2]], ["Glucose-6-Phosphate", "Glutamate"],
    tellheight = false, tellwidth = false, halign = :right, valign = :center, labelsize = 14)
fig

FileIO.save(joinpath(imgpath, "Glu_ED.pdf"), fig) 

CSV.write(joinpath(resultspath, "Results_ED.csv"), res[1])