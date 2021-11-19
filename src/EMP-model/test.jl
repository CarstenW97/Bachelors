using JuMP

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel

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

nadh_nad_ratio = 10

gln_model = GlnModel.gln_model(
    glc_ext = 50e-3,
    lac_ext = 1e-4,
    nh4_ext = 10e-3,
    ac_ext = 1e-16,
    etoh_ext = 1e-4,
    atp_adp_ratio = 10.0,
    nadh_nad_ratio = nadh_nad_ratio,
    num_ms = 10,
)

solve_1 = GlnModel.max_mu!(gln_model, prev_sol)

solve_2 = GlnModel.find_nearest!(gln_model, prev_sol) # find a solution close to the previous one, smoothen things

for sym in syms
    push!(eval(sym), value(gln_model[sym]))
end
# update old previous solution
GlnModel.update!(gln_model, prev_sol)
