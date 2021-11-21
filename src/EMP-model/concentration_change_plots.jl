using ColorSchemes, CairoMakie, FileIO, JuMP
using ProgressMeter

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel

# imgpath = joinpath("docs", "imgs", "EMP-model", "Concentration_changes")

#=
Warm up problem to get all the symbols in the model.
Create arrays using these symbols as the names.
This is the fast way of doing:
mu = Float64[]
atp = Float64[]
etc.
=#
test_model = GlnModel.gln_model();
syms = Symbol.(all_variables(test_model))
for sym in syms
    eval(:($sym=Float64[]))
end

prev_sol=Dict{Symbol, Float64}() # a variable to store the previous solution

con_lb = 1e-5
con_ub = 1e-7
@showprogress for var in range(log(con_lb); stop=log(con_ub), length=20)
    etoh_ext = exp(var)

    gln_model = GlnModel.gln_model(
        glc_ext = 50e-3,
        lac_ext = 1e-4,
        nh4_ext = 10e-3,
        ac_ext = 1e-16,
        etoh_ext = etoh_ext,
        # atp_adp_ratio = -1,
        # nadh_nad_ratio = -1, # unconstrained
        num_ms = 10,
    )

    solve_1 = GlnModel.max_mu!(gln_model, prev_sol)

    if !solve_1
        continue # if mu not found, skip
    end

    solve_2 = GlnModel.find_nearest!(gln_model, prev_sol) # find a solution close to the previous one, smoothen things

    if solve_1 && solve_2
        for sym in syms
            push!(eval(sym), value(gln_model[sym]))
        end
        # update old previous solution
        GlnModel.update!(gln_model, prev_sol)
    end
    
end

x = []
y = []
z = []

x = exp.(etoh_e)
y = exp.(atp)./exp.(adp)
z = exp.(nadh)./exp.(nad)

f = Figure()
ax = Axis(f[1,1], xscale=log10, 
          xlabel="External Ethanol concentration",
          ylabel="Relative ratio change",
)

scatter1 = scatter!(ax, x, y)
line1 = lines!(ax, x, y)
scatter2 = scatter!(ax, x, z)
line2 = lines!(ax, x, z)

f[1, 2] = Legend(
    f,
    [[scatter1, line1], [scatter2, line2]],
    ["ATP/ADP ratio", "NADH/NAD ratio"],
)

f

FileIO.save(joinpath(imgpath, "Ratio_change_Ethanol.pdf"), f)