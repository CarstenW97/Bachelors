using JuMP
using KNITRO

# Constraints (used the same for now)
UB = 100
LB = 1e-8

molenaar = Model(optimizer_with_attributes(KNITRO.Optimizer,
    "ms_enable" => 1,
    "opttol" => 1E-16,
    "opttolabs" => 1e-16,
    "honorbnds" => 1,
    "ms_maxsolves" => 2))
JuMP.set_silent(molenaar)

michaelis_menten(metabolite, protein, kcat, Km) =
        (kcat * protein * metabolite)/(Km + metabolite)
register(ancat, :michaelis_menten, 4, michaelis_menten; autodiff = true)

@variables molenaar begin


    # Concentrations
    LB <= sub_out <= UB # substarte outside concentration
    LB <= sub_in  <= UB # subtrate inside concentration
    LB <= ribo    <= UB # ribosome concentration
    LB <= lipid   <= UB # lipid concentration 
    LB <= prec    <= UB # precursor concentration
    LB <= transp  <= UB # transporter concentration
    LB <= metab   <= UB # metabolic enzyme concentration
    LB <= lip_syn <= UB # lipid biosynthesis enzyme concentration

    # Rates
    LB <= v_transp <= UB # transporter rate
    LB <= v_metab <= UB  # matabolic rate
    LB <= v_lip <= UB    # lipid biosynthesis rate
    LB <= v_ribo <= UB   # ribosomale rate

    # Cellular constants
    LB <= surface <= UB      # surface area of the cell 
    LB <= max_prot <= UB     # maximum protein concentration
    0 <= used_ribo[1:4] <= 1 # ribosomes used for synthesis
    LB <= l_to_t <= UB       # ratio of lipids to transporters
end

