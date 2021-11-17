include(joinpath("Thesis model", "Mos model.jl"))
import .Mos_Model

res_mo = []
for var in [0.01:0.05:5;]
    push!(res_mo, Mos_Model.mos_model(var))
end

################################################################

using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Thesis model", "Figures", "EMP_Model")
resultspath = joinpath("Thesis model", "Results", "EMP_Model")

res_df = Array{String,Float64}(res_emp)
df = DataFrame[]
for i = 1:100
    push!(df, res_emp[i])
end
CSV.write(joinpath(resultspath, "Results_EMP.csv"), res_emp)


x = [results_emp["mu"] for results_emp in res_emp]
a = [results_emp["g6p"] for results_emp in res_emp]
b = [results_emp["glu"] for results_emp in res_emp]

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
    halign = :right,
    valign = :center,
    labelsize = 14,
)
fig

FileIO.save(joinpath(imgpath, "Glu_EMP.pdf"), fig)
