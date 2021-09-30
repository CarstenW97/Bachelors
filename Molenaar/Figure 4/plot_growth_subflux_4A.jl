using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 4", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["sub_out"] for res in stats]
a = [res["cat_flux"] for res in stats]
b = [res["met_flux"] for res in stats]

fig = Figure()
ax = Axis(fig[1,1])
scatter1 = scatter!(ax, x, a)
line1 =lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 =lines!(ax, x, b)

ax.xlabel = "Substrate concentration"
ax.ylabel = "Fraction of the total substarte flux"

fig[1, 1] = Legend(fig, [[scatter1, line1], [scatter2, line2]], ["MetEf_flux", "CatEf_flux"],
    tellheight = false, tellwidth = false, halign = :right, valign = :center, labelsize = 14)
fig

FileIO.save(joinpath(imgpath, "Figure_4A.pdf"), fig)