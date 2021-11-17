using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 4", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["mu"] for res in stats]
a = [res["cat_flux"] for res in statsb]
b = [res["met_flux"] for res in statsb]

fig = Figure()
ax = Axis(fig[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)

ax.xlabel = "Growth rate"
ax.ylabel = "Fraction of the total substarte flux"
fig

FileIO.save(joinpath(imgpath, "Figure_4B.pdf"), fig)
