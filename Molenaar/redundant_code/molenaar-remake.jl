module MolenaarRemake
using JuMP
using KNITRO

function molenaar_remake(sub_out_con)
    
    molenaar = Model(optimizer_with_attributes(KNITRO.Optimizer,
        "ms_enable" => 1,
        "opttol" => 1E-16,
        "opttolabs" => 1e-16,
        "honorbnds" => 1,
        "ms_maxsolves" => 10))
    JuMP.set_silent(molenaar)

    # Constraints (used the same for now)
    UB = 100
    LB = 1e-8

    @variables molenaar begin
        LB <= trS   <= UB # substrate transporter
        LB <= rib   <= UB # ribosome
        LB <= catef <= UB # catalytic enzyme
        LB <= metef <= UB # metabolic enzyme
        LB <= prc   <= UB # precursor enzyme (activating enzyme)
        LB <= lpb   <= UB # lipid synthesis enzyme

        LB <= Si    <= UB # intracellular substrate
        LB <= M     <= UB # intermediate
        LB <= P     <= UB # precursor
        LB <= lip   <= UB # Lipid

        LB <= atp   <= UB
        LB <= adp   <= UB

        LB <= S     <= UB # substrate concentration outside
        
        LB <= mu <= UB

        0 <= alpha[1:6] <= 1 # fraction of ribosome engaged in synthesis

        LB <= beta <= UB

        LB <= v_trS    <= UB # transporter rate
        LB <= v_lpb    <= UB # lipid biosynthesis rate
        LB <= v_rib    <= UB # ribosomale rate
        LB <= v_catef  <= UB # catabolic enzyme rate
        LB <= v_metef  <= UB # metabolic enzyme rate
        LB <= v_prc    <= UB
    end

    @NLparameters molenaar begin
        kcat_trS   == 7
        kcat_prc   == 5
        kcat_rib   == 3
        kcat_catef == 10
        kcat_metef == 5
        kcat_lpb   == 5

        Km_trS   == 1
        Km_prc   == 1
        Km_rib   == 1
        Km_catef == 1
        Km_metef == 1
        Km_lpb   == 1

        Km_atp   == 1   # atp affinity of the activating enzyme prc
        Km_adp   == 0.5 # adp affinity of catef and metef
        gamm     == 0.6 # efficiency of catef 

        sA_trS == 1 # surface area of the transporter
        sA_lip == 1 # surface area of the lipids

        S == sub_out_con
    end

    @NLconstraints molenaar begin
        (v_trS - v_metef - v_catef) - mu * Si     == 0
        beta * ((sA_trS + sA_lip)*(lip + trS)) == 1

        sum(alpha[i] for i = 1:6)     == 1
        alpha[1] * v_rib - mu * trS   == 0
        alpha[2] * v_rib - mu * prc   == 0
        alpha[3] * v_rib - mu * catef == 0
        alpha[4] * v_rib - mu * metef == 0
        alpha[5] * v_rib - mu * rib   == 0
        alpha[6] * v_rib - mu * lpb   == 0

        (v_metef + gamm * v_catef - v_prc) - mu * atp == 0
        (v_prc - v_metef - gamm * v_catef) - mu * adp + mu * (atp + adp) == 0 

        v_rib   == (kcat_rib * rib * P) / (Km_rib + P)
        v_trS   == (kcat_trS * trS * S) / (Km_trS + S)
        v_prc   == (kcat_prc * prc * M * atp) / ((M * atp) + Km_prc * atp + Km_atp * M)
        v_catef == (kcat_catef * catef * Si * adp) / (Si * adp + Si * catef + adp * Km_catef)
        v_metef == (kcat_metef * metef * Si * adp) / (Si * adp + Si * metef + adp * Km_metef)
        v_lpb   == (kcat_lpb * lpb * P) / (Km_lpb + P)

        trS + rib + catef + metef + prc + lpb  == 1 # maximal intracellular protein concentration
        atp + adp                              == 1 # total concentration of energy intermediates

        (trS) <= lip * (lip/trS) # membrane integrity ((lip*trS) is the PLmax from the model)

    end

    @objective(molenaar, Max, mu)
    optimize!(molenaar)

    res = Dict(
        "S"       => value(S),
        "mu"      => value(mu),
        "beta"    => value(beta),
        "alpha"   => value.(alpha),
        "v_rib"   => value(v_rib),
        "v_trS"   => value(v_trS),
        "v_catef" => value(v_catef),
        "v_metef" => value(v_metef),
        "v_lpb"   => value(v_lpb),
        "atp"     => value(atp),
        "adp"     => value(adp))
    return res
end 

end #module