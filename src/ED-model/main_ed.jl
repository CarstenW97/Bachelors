include(joinpath("Bachelors", "src", "ED-model", "glu_ed_model.jl"))
import .Glu_ED_Model

res_ed = []
for var in [0.0001:0.0005:0.005;]
  push!(res_ed, Glu_ED_Model.glu_ed_model(var))
end
