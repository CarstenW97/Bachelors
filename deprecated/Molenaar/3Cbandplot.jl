using Colors, ColorSchemes
using CairoMakie, FileIO
using CSV, DataFrames
using Statistics

imgpath = joinpath("Molenaar", "Figure 3", "Plots")
resultspath = joinpath("Molenaar", "Results")

x = [res["mu"] for res in results]
y = [res["used_ribo"] for res in results]

a = [] # transporter rate
b = [] # matabolic rate
c = [] # ribosomale rate
d = [] # lipid biosynthesis rate

for i = 1:100
    push!(a, y[i][1])
    push!(b, y[i][2])
    push!(c, y[i][3])
    push!(d, y[i][4])
end

a = Float64.(a)
b = Float64.(b)
c = Float64.(c)
d = Float64.(d)

metabolism = b
ribosomes = c
transporters = a
lipids = d

fracs = [metabolism, ribosomes, transporters, lipids]
fraclabels = ["Metabolism", "Ribosomes", "Transporters", "Lipidsynthesis"]
function plot_proteome(ax, x, fracs, fraclabels)
    d = cumsum(hcat(fracs...)'; dims = 1)
    d = [zeros(size(d, 2))'; d ./ d[end, :]']
    for i = 1:length(fraclabels)
        band!(ax, x, d[i, :], d[i+1, :], label = fraclabels[i])
    end
end

f = Figure();
ax = Axis(f[1, 1]);
plot_proteome(ax, 1:100, fracs, fraclabels)
ax.xlabel = "Simulations with rising substarte concentration"
ax.ylabel = "Proteome fraction"
f[1, 2] = Legend(f, ax, "Protein", unique = true, framevisible = false)
f

FileIO.save(joinpath(imgpath, "Figure_3Cband.pdf"), f)
