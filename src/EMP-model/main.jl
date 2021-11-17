using JuMP
using JSON

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel

#=
Solve the optimization problem
=#
varname = "lac_ext"

!isdir(joinpath("src", "EMP-model", varname)) && mkdir(joinpath("src", "EMP-model", varname))

for var in range(1e-3, 1e-1; length=10)

    m = GlnModel.gln_model(
        glc_ext = 0.05,
        lac_ext = var,
        nh4_ext = 0.01,
        ac_ext = 1e-16,
        etoh_ext = 1e-6,
        atp_adp_ratio = 10.0,
        nadh_nad_ratio = 0.2,
    )

    isnothing(m) && continue

    open(joinpath("src", "EMP-model", varname, "sol_with_$(varname)_set_$(string(var)).json"), "w") do io
        JSON.print(io, Dict(string(id) => value(m[id]) for id in Symbol.(all_variables(m))))
    end

end
