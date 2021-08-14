using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["mu"] for res in results]
y = [res["used_ribo"] for res in results]

a = [] # transporter rate
b = [] # matabolic rate
c = [] # ribosomale rate
d = [] # lipid biosynthesis rate

for i in 1:20
    push!(a, y[i][1])
    push!(b, y[i][2])
    push!(c, y[i][3])
    push!(d, y[i][4])
end

a = Float64.(a)
b = Float64.(b)
c = Float64.(c)
d = Float64.(d)

fig = Figure()
ax = Axis(fig[1,1])
scatter1 = scatter!(ax, x, a) # blue
line1 = lines!(ax, x, a)
scatter2 = scatter!(ax, x, b) # yellow
line2 = lines!(ax, x, b)
scatter3 = scatter!(ax, x, c) # green
line3 = lines!(ax, x, c)
scatter4 = scatter!(ax, x, d) # pink
line4 = lines!(ax, x, d)

ax.xlabel = "Growth rate"
ax.ylabel = "Protein fraction"
fig[1, 1] = Legend(fig, [[scatter1, line1], [scatter2, line2], [scatter3, line3], [scatter4, line4]], ["Transporter", "Metabolism", "Ribosomes", "Lipids"],
    tellheight = false, tellwidth = false, halign = :right, valign = :top, labelsize = 14)
fig

FileIO.save(joinpath(imgpath, "Figure_3C.pdf"), fig)
  