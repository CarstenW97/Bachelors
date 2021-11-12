using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("docs", "imgs", "EMP-model")
resultspath = joinpath("docs", "results", "EMP-model")

x = [exp(results_emp["glc_ext"]) for results_emp in res_emp]
y = [exp(results_emp["pts"]) for results_emp in res_emp]

fig = Figure()
ax = Axis(fig[1,1])
scatter_glc = scatter!(ax, x, y)
line_glc = lines!(ax, x, y)

ax.xlabel = "External glucose concetration [mM]"
ax.ylabel = "PTS concentration [g enz / g DW]"

fig

FileIO.save(joinpath(imgpath, "EMP_PTScon_Glc.pdf"), fig)
