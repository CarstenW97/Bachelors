include(joinpath("src", "EMP-model", "glu_emp_model.jl"))
import .Glu_EMP_Model

res_emp = []
for var in [0.0001:0.0005:0.005;]
  push!(res_emp, Glu_EMP_Model.glu_emp_model(var))
end

#Glu_EMP_Model.glu_emp_model(0.05)
