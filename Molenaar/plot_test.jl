using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["sub_out"] for res in results_test]
y = [res["mu"] for res in results_test]

fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax, x, y)
lines!(ax, x, y)

ax.xlabel = "Substrate concentration outside"
ax.ylabel = "Growth rate"
fig