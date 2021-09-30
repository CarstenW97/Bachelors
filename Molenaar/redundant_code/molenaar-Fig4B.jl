module MolenaarAltB
using JuMP
using KNITRO

function molenaar_altb(sub_out_con)
    
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

    michaelis_menten_adp(sub1, sub2, prot1, prot2, kcat, Km1, Km2) = # sub1=sub_in; sub2=adp; Km1=protein; Km2=adp, prot1=catef, prot2=metef
            (kcat * sub1 * sub2 * prot1 * prot2) / ((sub1 * sub2) + (sub1 * Km2)+ (sub2 * Km1))
    register(molenaar, :michaelis_menten_adp, 7, michaelis_menten_adp; autodiff = true)

    michaelis_menten_atp(sub1, sub2, protein, kcat, Km1, Km2) = # sub1=inter; sub2=atp; Km1=activ; Km2=atp
            (kcat *sub1 * sub2 * protein) / ((sub1 * sub2) + (sub1 * Km2)+ (sub2 * Km1))
    register(molenaar, :michaelis_menten_atp, 6, michaelis_menten_atp; autodiff = true)

    # Constraints (used the same for now)
    UB = 100
    LB = 1e-8

    @variables molenaar begin
        # growth rate
        LB <= mu <= UB

        # Concentrations
        LB <= sub_out <= UB # substarte outside concentration
        LB <= sub_in  <= UB # subtrate inside concentration
        LB <= lipid   <= UB # lipid concentration 
        LB <= prec    <= UB # precursor concentration

        LB <= inter   <= UB # intermediate concentration
        LB <= atp     <= UB # atp concentration
        LB <= adp     <= UB # adp concentration

        LB <= ribo    <= UB # ribosome concentration
        LB <= transp  <= UB # transporter concentration
        LB <= lip_syn <= UB # lipid biosynthesis enzyme concentration
        LB <= catef   <= UB # catabolic ineffecient enzyme
        LB <= metef   <= UB # metabolic effecient enzyme

        LB <= activ   <= UB # activating enzyme

        # Rates
        LB <= v_transp <= UB # transporter rate
        LB <= v_lip    <= UB # lipid biosynthesis rate
        LB <= v_ribo   <= UB # ribosomale rate
        LB <= v_catef  <= UB # catabolic enzyme rate
        LB <= v_metef  <= UB # metabolic enzyme rate

        LB <= v_activ  <= UB # activating enzyme rate

        # Cellular constants
        LB <= beta           <= UB # surface to area ratio  
        0  <= used_ribo[1:6] <= 1  # ribosomes used for synthesis
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
        kcat_catef  == 10
        Km_catef    == 1
        kcat_metef  == 5
        Km_metef    == 1
        
        kcat_activ  == 1
        Km_activ    == 1
        gamma       == 0.6 # amount of atp created by catef per reaction
        Km_atp      == 1
        Km_adp      == 0.5
    end

    @NLconstraints molenaar begin
        v_ribo   == michaelis_menten(prec, ribo, kcat_ribo, Km_ribo)
        v_transp == michaelis_menten(sub_out, transp, kcat_transp, Km_transp)
        v_lip    == michaelis_menten(prec, lip_syn, kcat_lip, Km_lip)
        v_catef  == michaelis_menten_adp(sub_in, adp, catef, metef, kcat_catef, Km_catef, Km_adp) 
        v_metef  == michaelis_menten_adp(sub_in, adp, catef, metef, kcat_metef, Km_metef, Km_adp) 
        v_activ  == michaelis_menten_atp(inter, atp, activ, kcat_activ, Km_activ, Km_atp)  

        beta * (lipid + transp)                == 1     # membrane proportions
        sum(used_ribo[i] for i=1:6)            == 1     # ribosome proportions
        catef + metef + ribo + lip_syn + activ <= 1.0   # intracellular density
        transp                                 <= lipid # transporter density

        atp + adp                              == 1

        # ribosome activity producing enzymes
        used_ribo[1] * v_ribo - mu * transp  == 0
        used_ribo[2] * v_ribo - mu * ribo    == 0
        used_ribo[3] * v_ribo - mu * lip_syn == 0
        used_ribo[4] * v_ribo - mu * catef   == 0
        used_ribo[5] * v_ribo - mu * metef   == 0
        used_ribo[6] * v_ribo - mu * activ   == 0

        # metabolic enzyme balances
        (v_transp - v_metef - v_catef) - mu * sub_in      == 0
        (v_metef + v_catef - v_activ) - mu * inter        == 0  
        (v_metef + gamma * v_catef - v_activ) - mu * atp  == 0  
        (v_activ - v_lip - v_ribo) - mu * prec            == 0  
        (v_activ - v_metef - gamma * v_catef) - mu * adp  == 0  
        v_lip - mu * lipid                                == 0
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