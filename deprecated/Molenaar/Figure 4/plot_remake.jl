using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 4", "Plots")
resultspath = joinpath("Molenaar", "Results")

#x = [res["S"] for res in results_test]
#a = [res["v_catef"] for res in results_test]
#b = [res["v_metef"] for res in results_test]

#fig = Figure()
#ax = Axis(fig[1,1])
#scatter1 = scatter!(ax, x, a)
#line1 =lines!(ax, x, a)
#scatter2 = scatter!(ax, x, b)
#line2 =lines!(ax, x, b)

#ax.xlabel = "Growth rate"
#ax.ylabel = "Fraction of the total substarte flux"
#fig

#FileIO.save(joinpath(imgpath, "Figure_remake.pdf"), fig)

# Remake 2

x = [res["sub_out"] for res in results_test]
a = [res["v_catef"] for res in results_test]
b = [res["v_metef"] for res in results_test]

fig = Figure()
ax = Axis(fig[1, 1])
scatter1 = scatter!(ax, x, a)
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b)
line2 = lines!(ax, x, b)

ax.xlabel = "Substrate concentration"
ax.ylabel = "Fraction of the total substarte flux"

fig[1, 1] = Legend(
    fig,
    [[scatter1, line1], [scatter2, line2]],
    ["MetEf_flux", "CatEf_flux"],
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :center,
    labelsize = 14,
)
fig

FileIO.save(joinpath(imgpath, "Figure_remake2.pdf"), fig)
