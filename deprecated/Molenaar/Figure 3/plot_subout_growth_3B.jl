using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 3", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["sub_out"] for res in results]
y = [res["mu"] for res in results]

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, x, y)
lines!(ax, x, y)

ax.xlabel = "Substrate concentration outside"
ax.ylabel = "Growth rate (1/h)"
fig

FileIO.save(joinpath(imgpath, "Figure_3B.pdf"), fig)
