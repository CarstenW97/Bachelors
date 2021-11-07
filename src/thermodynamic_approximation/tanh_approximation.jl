using CairoMakie, LsqFit

#=
Use tanh to approximate the thermodynamic constraint
flux_forward = kinetics * (1 - exp(ΔG/RT))
flux_backward = kinetics * (exp(-ΔG/RT) - 1)
=#

xs = [
    range(-10, stop = -3, length = 1000)
    range(-3, stop = 3, length = 10000)
    range(3, stop = 10, length = 1000)
];

ys = [x <= 0 ? 1.0 - exp(x) : exp(-x) - 1.0 for x in xs]

@. fitfunc(x, p) = tanh(p[1] * x)
p0 = [1.0]

fit = curve_fit(fitfunc, xs, ys, p0)
vs = round.(coef(fit), digits = 1)

ysfit = fitfunc.(xs, Ref(vs));

fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "ΔG/RT", ylabel = "Thermodynamic factor")
lines!(ax1, xs, ys, label = "True")
lines!(ax1, xs, ysfit, label = "Approximation")
axislegend()

ax2 = Axis(fig[1, 2], xlabel = "ΔG/RT")
lines!(ax2, xs, ys)
lines!(ax2, xs, ysfit)
xlims!(ax2, -0.3, 0.3)
ylims!(ax2, -0.3, 0.3)

ax3 = Axis(fig[2, 1:2], xlabel = "ΔG/RT", ylabel = "% Error")
lines!(ax3, xs, (ys .- ysfit) ./ ys .* 100)

fig
CairoMakie.FileIO.save(joinpath("docs", "imgs", "tanh_approximation.pdf"), fig)
