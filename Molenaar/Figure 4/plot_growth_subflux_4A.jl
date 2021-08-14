using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 4", "Plots")
resultspath = joinpath("Molenaar", "Results")

fig = Figure()
ax = Axis(fig[1,1])
scatter1 = scatter!(ax, x, y)
line1 =lines!(ax, x, y)
scatter2 = scatter!(ax, x, z)
line2 =lines!(ax, x, z)


ax.xlabel = ""
ax.ylabel = ""
fig


FileIO.save(joinpath(imgpath, "Figure_4A.pdf"), fig)