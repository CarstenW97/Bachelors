using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 3", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["mu"] for res in results]
y = [res["v_ribo"] for res in results]

fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax, x, y)
lines!(ax, x, y)

ax.xlabel = "Growth rate"
ax.ylabel = "Ribosome synthesis rate"
fig

FileIO.save(joinpath(imgpath, "Figure_3E.pdf"), fig)