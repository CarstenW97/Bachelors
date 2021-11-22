using JuMP

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel


prev_sol=Dict{Symbol, Float64}() # a variable to store the previous solution

nadh_nad_ratio = 0.1

gln_model = GlnModel.gln_model(
    glc_ext = 50e-3,
    lac_ext = 1e-4,
    nh4_ext = 10e-3,
    ac_ext = 1e-16,
    etoh_ext = 1e-4,
    atp_adp_ratio = 10,
    nadh_nad_ratio = nadh_nad_ratio,
    num_ms = 10,
    silence=false,
)

solve_1 = GlnModel.max_mu!(gln_model, prev_sol)

solve_2 = GlnModel.find_nearest!(gln_model, prev_sol) # find a solution close to the previous one, smoothen things


# update old previous solution
GlnModel.update!(gln_model, prev_sol)
