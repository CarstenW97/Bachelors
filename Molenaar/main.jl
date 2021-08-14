include(joinpath("Molenaar", "molenaar-basic.jl"))
import .MolenaarBasic

res = MolenaarBasic.molenaar_basic(2)


# Figures
results = []
for sext in 0.001:20
  push!(results, MolenaarBasic.molenaar_basic(sext))
end



# Test for more data points between 1 and 2
j=0.1
results_test = []
for sext in 0.001:20:j
  push!(results_test, MolenaarBasic.molenaar_basic(sext))
  if sext == 2
    j = 1
  end
end