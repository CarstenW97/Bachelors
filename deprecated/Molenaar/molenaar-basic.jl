module MolenaarBasic
using JuMP
using KNITRO

function molenaar_basic(sub_out_con)
    molenaar = Model(
        optimizer_with_attributes(
            KNITRO.Optimizer,
            "ms_enable" => 1,
            "opttol" => 1E-16,
            "opttolabs" => 1e-16,
            "honorbnds" => 1,
            "ms_maxsolves" => 10,
        ),
    )
    JuMP.set_silent(molenaar)

    michaelis_menten(metabolite, protein, kcat, Km) =
        (kcat * protein * metabolite) / (Km + metabolite)
    register(molenaar, :michaelis_menten, 4, michaelis_menten; autodiff = true)

    # Constraints (used the same for now)
    UB = 100
    LB = 1e-8

    @variables molenaar begin
        # growth rate
        LB <= mu <= UB

        # Concentrations
        LB <= sub_out <= UB # substarte outside concentration
        LB <= sub_in <= UB # subtrate inside concentration
        LB <= ribo <= UB # ribosome concentration
        LB <= lipid <= UB # lipid concentration 
        LB <= prec <= UB # precursor concentration
        LB <= transp <= UB # transporter concentration
        LB <= metab <= UB # metabolic enzyme concentration
        LB <= lip_syn <= UB # lipid biosynthesis enzyme concentration

        # Rates
        LB <= v_transp <= UB # transporter rate
        LB <= v_metab <= UB  # matabolic rate
        LB <= v_lip <= UB    # lipid biosynthesis rate
        LB <= v_ribo <= UB   # ribosomale rate

        # Cellular constants
        LB <= beta <= UB         # surface to area ratio  
        0 <= used_ribo[1:4] <= 1 # ribosomes used for synthesis
    end

    @NLparameters molenaar begin
        kcat_ribo == 3
        Km_ribo == 1
        kcat_transp == 7
        Km_transp == 1
        kcat_metab == 5
        Km_metab == 1
        kcat_lip == 5
        Km_lip == 1
        sub_out == sub_out_con
    end

    @NLconstraints molenaar begin
        v_ribo == michaelis_menten(prec, ribo, kcat_ribo, Km_ribo)
        v_transp == michaelis_menten(sub_out, transp, kcat_transp, Km_transp)
        v_metab == michaelis_menten(sub_in, metab, kcat_metab, Km_metab)
        v_lip == michaelis_menten(prec, lip_syn, kcat_lip, Km_lip)

        beta * (lipid + transp) == 1     # membrane proportions
        sum(used_ribo[i] for i = 1:4) == 1   # ribosome proportions
        metab + ribo + lip_syn <= 1.0      # intracellular density
        transp <= lipid                    # transporter density

        # ribosome activity producing enzymes
        used_ribo[1] * v_ribo - mu * transp == 0
        used_ribo[2] * v_ribo - mu * metab == 0
        used_ribo[3] * v_ribo - mu * ribo == 0
        used_ribo[4] * v_ribo - mu * lip_syn == 0

        # metabolic enzyme balances
        (v_transp - v_metab) - mu * sub_in == 0
        (v_metab - v_lip - v_ribo) - mu * prec == 0
        v_lip - mu * lipid == 0
    end

    @objective(molenaar, Max, mu)
    optimize!(molenaar)

    res = Dict(
        "mu" => value(mu),
        "used_ribo" => value.(used_ribo),
        "sub_out" => value(sub_out),
        "sub_in" => value(sub_in),
        "v_ribo" => value(v_ribo),
    )
    return res
end

end #module
