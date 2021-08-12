using JuMP
using KNITRO

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

# Constraints (used the same for now)
UB = 100
LB = 1e-8

@variables molenaar begin
    # growth rate
    LB <= mu <= UB

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

@NLparamater molenaar begin
    kcat_ribo == 3
    Km_ribo == 1
end

@NLconstraints molenaar begin
    v_ribo == michaelis_menten(prec, ribo, kcat_ribo, Km_ribo)
    surface * (lip_syn + transp) == 1                           # membrane proportions
    sum(alpha[i] for i=1:4) == 1                                # ribosome proportions
    metab + ribo + lip_syn <= 1.0                               # intracellular density
    ts + tw <= L                                                # transporter density

    # ribosome activity producing enzymes
    alpha[1] * v_ribo - mu * transp == 0
    alpha[2] * v_ribo - mu * metab == 0
    alpha[3] * v_ribo - mu * ribo == 0
    alpha[4] * v_ribo - mu * lip_syn == 0

    # metabolic enzyme balances
    (v_ts - v_c - v_a) - mu * sub_in == 0
    (v_a - v_l - v_r) - mu * prec == 0
    v_l - mu * lip_syn == 0
end