using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 3", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["mu"] for res in results]
y = [res["sub_in"] for res in results]

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, x, y)
lines!(ax, x, y)

ax.xlabel = "Growth rate (1/h)"
ax.ylabel = "Substrate concentration inside"
fig

FileIO.save(joinpath(imgpath, "Figure_3D.pdf"), fig)
