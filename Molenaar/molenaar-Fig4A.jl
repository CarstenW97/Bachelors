module MolenaarAlt
using JuMP
using KNITRO

function molenaar_alt(sub_out_con)
    
    molenaar = Model(optimizer_with_attributes(KNITRO.Optimizer,
        "ms_enable" => 1,
        "opttol" => 1E-16,
        "opttolabs" => 1e-16,
        "honorbnds" => 1,
        "ms_maxsolves" => 10))
    JuMP.set_silent(molenaar)

    michaelis_menten(metabolite, protein, kcat, Km) =
            (kcat * protein * metabolite)/(Km + metabolite)
    register(molenaar, :michaelis_menten, 4, michaelis_menten; autodiff = true)

    # Constraints (used the same for now)
    UB = 100
    LB = 1e-8

    @variables molenaar begin
        # growth rate
        LB <= mu <= UB

        # Concentrations
        LB <= sub_out <= UB # substrate outside concentration
        LB <= sub_in  <= UB # subtrate inside concentration
        LB <= ribo    <= UB # ribosome concentration
        LB <= lipid   <= UB # lipid concentration 
        LB <= prec    <= UB # precursor concentration
        LB <= transp  <= UB # transporter concentration
        LB <= lip_syn <= UB # lipid biosynthesis enzyme concentration

            # new
        LB <= catef   <= UB # catabolic ineffecient enzyme
        LB <= metef   <= UB # metabolic effecient enzyme

        # Rates
        LB <= v_transp <= UB # transporter rate
        LB <= v_lip    <= UB # lipid biosynthesis rate
        LB <= v_ribo   <= UB # ribosomale rate

            # new
        LB <= v_catef  <= UB # catabolic enzyme rate
        LB <= v_metef  <= UB # metabolic enzyme rate

        # Cellular constants
        LB <= beta           <= UB # surface to area ratio  
        0  <= used_ribo[1:5] <= 1  # ribosomes used for synthesis
    end

    @NLparameters molenaar begin
        kcat_ribo   == 3
        Km_ribo     == 1
        kcat_transp == 7
        Km_transp   == 1
        kcat_metab  == 5
        Km_metab    == 1
        kcat_lip    == 5
        Km_lip      == 1
        sub_out     == sub_out_con

            #new
        kcat_catef  == 10
        Km_catef    == 1
        kcat_metef  == 5
        Km_metef    == 1
        
        gamma       == 0.8
    end

    @NLconstraints molenaar begin
        v_ribo == michaelis_menten(prec, ribo, kcat_ribo, Km_ribo)
        v_transp == michaelis_menten(sub_out, transp, kcat_transp, Km_transp)
        v_lip == michaelis_menten(prec, lip_syn, kcat_lip, Km_lip)

            #new
        v_catef == michaelis_menten(sub_in, catef, kcat_catef, Km_catef)
        v_metef == michaelis_menten(sub_in, metef, kcat_metef, Km_metef)

        beta * (lipid + transp)        == 1     # membrane proportions
        sum(used_ribo[i] for i=1:5)    == 1     # ribosome proportions
        catef + metef + ribo + lip_syn <= 1.0   # intracellular density
        transp                         <= lipid # transporter density

        # ribosome activity producing enzymes
        used_ribo[1] * v_ribo - mu * transp  == 0
        used_ribo[2] * v_ribo - mu * ribo    == 0
        used_ribo[3] * v_ribo - mu * lip_syn == 0
        used_ribo[4] * v_ribo - mu * catef   == 0
        used_ribo[5] * v_ribo - mu * metef   == 0

        # metabolic enzyme balances
        (v_transp - v_metef - v_catef) - mu * sub_in     == 0
        (v_metef + gamma * v_catef - v_lip - v_ribo) - mu * prec == 0
        v_lip - mu * lipid                     == 0
    end

    @objective(molenaar, Max, mu)
    optimize!(molenaar)

    res = Dict(
            "mu"       => value(mu),
            "met_flux" => value(v_metef),
            "cat_flux" => value(v_catef),
        )
    return res
end 

end #module
