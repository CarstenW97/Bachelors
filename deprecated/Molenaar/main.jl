include(joinpath("Molenaar", "molenaar-basic.jl"))
import .MolenaarBasic

res = MolenaarBasic.molenaar_basic(2)


# Figures
results = []
for sext in [0.01:0.05:5;]
  push!(results, MolenaarBasic.molenaar_basic(sext))
end

# Test for the remake of the molenaar model
include(joinpath("Molenaar", "molenaar_remake2.jl"))
import .MolenaarRemake2

#res = MolenaarRemake2.molenaar_remake2(5.0)

results_test = []
for sext in [0.01:0.05:5;]
  push!(results_test, MolenaarRemake2.molenaar_remake2(sext))
end

#########

include(joinpath("Molenaar", "molenaar-Fig4A.jl"))
import .MolenaarAlt

res = MolenaarAlt.molenaar_alt(2)

# Figure 4A
stats = []
for sext in [0.01:0.05:5;]
  push!(stats, MolenaarAlt.molenaar_alt(sext))
end

#########

include(joinpath("Molenaar", "molenaar-Fig4B.jl"))
import .MolenaarAltB

res = MolenaarAltB.molenaar_altb(2)

# Figure 4B
statsb = []
for sext in 0.001:20
  push!(statsb, MolenaarAltB.molenaar_altb(sext))
end