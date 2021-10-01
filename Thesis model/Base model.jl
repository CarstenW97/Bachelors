module BaseModel
using JuMP
using KNITRO

S_ext = 

function basemodel_1(S_ext)
    
    basemodel = Model(optimizer_with_attributes(KNITRO.Optimizer,
        "ms_enable" => 1,
        "opttol" => 1E-16,
        "opttolabs" => 1e-16,
        "honorbnds" => 1,
        "ms_maxsolves" => 10))
    JuMP.set_silent(basemodel)

    # Constraints (used the same for now)
    UB = 100
    LB = 1e-8

    @variables molenaar begin
        LB <= rib       <= UB # ribosomes
        LB <= trS       <= UB # substrate transporters
        LB <= EMP_enz   <= UB # EMP pathway enzyme
        LB <= ED_enz    <= UB # ED pathway enzyme
        LB <= lpd       <= UB # Lipidsynthesis enzyme 

        #LB <= S_ext     <= UB # extracellular substrate
        LB <= S_int     <= UB # intracellular substrate
        LB <= W         <= UB # EMP/ED product
        LB <= Lip       <= UB # Lipid

        LB <= atp       <= UB
        LB <= adp       <= UB
        LB <= nad       <= UB
        LB <= nadh      <= UB
 
        LB <= mu        <= UB # growth rate
        0 <= alpha[1:5] <= 1 # fraction of ribosome engaged in synthesis
        LB <= beta      <= UB # surface to area ratio

        LB <= v_rib     <= UB # ribosomale rate
        LB <= v_trS     <= UB # transporter rate
        LB <= v_EMP     <= UB # EMP enzyme rate
        LB <= v_ED      <= UB # ED enzyme rate
        LB <= v_lpb     <= UB # lipid biosynthesis rate
    end

    @NLparameters molenaar begin
        
        R = 8.3145 # gas constants
        T = 298.15 # tempreture [K]

        # enzyme rates
        kcat_rib   ==
        kcat_trS   == 
        kcat_EMP   == 
        kcat_ED    == 
        kcat_lpb   == 

        # Michaelis-menten constants
        #Km_rib     == 
        #Km_trS     == 
        #Km_EMP     == 
        #Km_ED      == 
        #Km_lpd     == 
        #Km_Sext    ==
        #Km_Sint    ==
        #Km_W       ==

        #Km_atp   ==    # atp affinity of the anabolic enzymes
        #Km_adp   ==    # adp affinity of EMP/ED enzyme
        #gamm     ==    # efficiency of ED

        sA_trS == 1 # surface area of the transporter
        sA_lip == 1 # surface area of the lipids

        # Gibbs energy of reaction
        dG_rib  ==
        dG_trS  ==
        dG_EMP  ==
        dG_ED   ==
        dG_lip  == 
        dG_Sint ==
        dG_W    == 0
    end

    @NLconstraints molenaar begin
        (v_trS - v_EMP - v_ED) - mu * S_int     == 0
        (v_EMP + v_ED - v_rib - v_lpb) - mu * W == 0
        (v_lpb) - mu * lip                      == 0

        beta * (lip + trS) == 1

        sum(alpha[i] for i = 1:6)       == 1
        alpha[1] * v_rib - mu * rib     == 0
        alpha[2] * v_rib - mu * trS     == 0
        alpha[3] * v_rib - mu * enz_EMP == 0
        alpha[4] * v_rib - mu * enz_ED  == 0
        alpha[5] * v_rib - mu * lpd     == 0

        (v_metef + gamm * v_catef - v_prc) - mu * atp == 0
        (v_prc - v_metef - gamm * v_catef) - mu * adp + mu * (atp + adp) == 0

        # v = [E] * kcat * (1-e^dG/R*T)
        v_rib == rib * kcat_rib * (1-e^(dG_rib/(R*T)))
        v_trS == trS * kcat_trS * (1-e^(dG_Sint/(R*T)))
        v_EMP == EMP_enz * kcat_EMP * (1-e^(dG_EMP/(R*T)))
        v_ED  == ED_enz * kcat_ED * (1-e^(dG_ED/(R*T)))
        v_lpd == lpd * kcat_lpd * (1-e^(dG_lip/(R*T)))

        rib + trS + EMP + ED + lpd       == 1 # maximal intracellular protein concentration
        atp + adp                        == 1 # total concentration of energy intermediates

        trS <= lip  # membrane integrity ((lip*trS) is the PLmax from the model)

    end

    @objective(basemodel, Max, mu)
    optimize!(basemodel)
    value(mu)

    res = Dict(
        "S_ext"  => value(S_ext),
        "S_int"  => value(S_int),
        "mu"     => value(mu),
        "beta"   => value(beta),
        "alpha"  => value.(alpha),
        "v_rib"  => value(v_rib),
        "v_trS"  => value(v_trS),
        "v_EMP"  => value(v_EMP),
        "v_ED"   => value(v_ED),
        "v_lpb"  => value(v_lpb),
        "atp"    => value(atp),
        "adp"    => value(adp))
    return res
end     

end #module