# module Glu_ED_Model
using JuMP
using KNITRO

# function glu_ed_model(glc_ext_input)
    glc_ext_input = 0.05
    glu_ed = Model(optimizer_with_attributes(KNITRO.Optimizer,
        "ms_enable" => 1,
        "ms_maxsolves" => 10))
    # JuMP.set_silent(glu_ed)

    # Thermodynamic function
    thermo_factor(x) = tanh(-0.7 * x) # x = ΔrG/RT
    register(glu_ed, :thermo_factor, 1, thermo_factor; autodiff = true)

    # Bounds
    LB_enz =    0.0     # [g enz / g DW cell]
    UB_enz =    1.0     # [g enz / g DW cell]
    LB_met = log(1e-6)  # [log(1 uM)]
    UB_met = log(10e-3) # [log(10 mM)]
    LB_v   = - 50       # [mmol/gDW/h]
    UB_v   =   50       # [mmol/gDW/h]
    LB_dg  = -100       # [kJ/mol]
    UB_dg  =  100       # [kJ/mol]

    @variables glu_ed begin
        # biomass [1/h] (g6p + ATP -> ADP)
        0 <= mu <= 10

        # Enzymes concentrations [g enz / g DW]
        LB_enz <= pts     <= UB_enz # enzyme for: glc external + ATP -> ADP + g6p
        LB_enz <= ed      <= UB_enz # enzyme for: g6p + NAD + ADP + NADP -> NADH + ATP + NADPH + pep + pyr
        LB_enz <= pyk     <= UB_enz # enzyme for: pep + ADP -> ATP + pyr
        LB_enz <= ldh     <= UB_enz # enzyme for: pyr + NADH -> NAD + lac
        LB_enz <= ppc     <= UB_enz # enzyme for: pep + CO2 -> oac + P
        LB_enz <= akgsyn  <= UB_enz # enzyme for: pyr + oac -> akg
        LB_enz <= gdhm    <= UB_enz # enzyme for: akg + NH4 + NADH -> glu + NAD
        LB_enz <= burn    <= UB_enz # enzyme for: ATP -> ADP
        LB_enz <= nadtrdh <= UB_enz # enzyme for: NAD + NADPH -> NADH + NADP
        LB_enz <= lp      <= UB_enz # enzyme for: lac internal -> lac external

        # Metabolite concentrations [log(M)]
        LB_met <= g6p   <= UB_met # glucose-6-phosphate
        LB_met <= pyr   <= UB_met # pyruvate
        LB_met <= pep   <= UB_met # phosphoenolpyrovate
        LB_met <= oac   <= UB_met # oxalacetate
        LB_met <= akg   <= UB_met # alpha-ketogluterate
        LB_met <= glu   <= UB_met # glutamate
        LB_met <= lac   <= UB_met # lactate
        LB_met <= atp   <= UB_met # ATP
        LB_met <= adp   <= UB_met # ADP
        LB_met <= nad   <= UB_met # NAD
        LB_met <= nadh  <= UB_met # NADH
        LB_met <= nadp  <= UB_met # NADP
        LB_met <= nadph <= UB_met # NADPH
        LB_met <= co2   <= UB_met
        LB_met <= phos  <= UB_met # phosphate
        LB_met <= nh4   <= UB_met

        # Fluxes [mmol/gDW/h]
        LB_v <= v_pts     <= UB_v
        LB_v <= v_ed      <= UB_v
        LB_v <= v_pyk     <= UB_v
        LB_v <= v_ldh     <= UB_v
        LB_v <= v_ppc     <= UB_v
        LB_v <= v_akgsyn  <= UB_v
        LB_v <= v_gdhm    <= UB_v
        LB_v <= v_burn    <= UB_v
        LB_v <= v_nadtrdh <= UB_v
        LB_v <= v_lp      <= UB_v

        # Thermodynamic variables
        LB_dg <= dg_pts     <= UB_dg
        LB_dg <= dg_ed      <= UB_dg
        LB_dg <= dg_pyk     <= UB_dg
        LB_dg <= dg_ldh     <= UB_dg
        LB_dg <= dg_ppc     <= UB_dg
        LB_dg <= dg_akgsyn  <= UB_dg
        LB_dg <= dg_gdhm    <= UB_dg
        LB_dg <= dg_burn    <= UB_dg
        LB_dg <= dg_nadtrdh <= UB_dg
        LB_dg <= dg_lp      <= UB_dg
    end

    @NLparameters glu_ed begin

        #R  == 8.3145e-3 # gas constants [kJ/K/mol]
        #T  == 298.15    # tempreature   [K]
        #RT == R * T
        RT == 8.3145e-3 * 298.15

        # Enzyme rates
        kcat_pts     == 213.75
        kcat_ed      == 268
        kcat_pyk     ==  58.4
        kcat_ldh     ==  31
        kcat_ppc     == 540
        kcat_akgsyn  ==   4
        kcat_gdhm    ==  20
        kcat_burn    ==  10
        kcat_nadtrdh == 167.9
        kcat_lp      == 100

        # Gibbs energy of reaction
        dG0_pts     == -16.7
        dG0_ed      == -81.8
        dG0_pyk     == -31.7
        dG0_ldh     == -23.7
        dG0_ppc     == -40.3
        dG0_akgsyn  == -60.7
        dG0_gdhm    == -33.4
        dG0_burn    == -57
        dG0_nadtrdh ==   0
        dG0_lp      == -10.3

        # media conditions
        glc_ext == log(glc_ext_input) # [log(50 mM)]
        lac_ext == log(10e-6) # [log(10 uM)]

        # capacity constraint
        total_proteome_mass_fraction == 0.26 # [g enz/gDW]

        # maintenance requirement
        min_burn_flux == 0.01 # [mmol/gDW/h]
    end

    @NLconstraints glu_ed begin
        # ΔG at concentrations
        dg_pts     == dG0_pts + RT * (adp + g6p - atp - glc_ext)
        dg_ed      == dG0_ed + RT * (nadph + nadh + atp + pep + pyr - adp -  nadp - nad - g6p)
        dg_pyk     == dG0_pyk + RT * (pyr + atp - adp - pep)
        dg_ldh     == dG0_ldh + RT * (lac + nad - nadh - pyr)
        dg_ppc     == dG0_ppc + RT * (pep - oac)
        dg_akgsyn  == dG0_akgsyn  + RT * (pyr + oac - akg)
        dg_gdhm    == dG0_gdhm + RT * (akg + nadh - glu - nad)
        dg_burn    == dG0_burn + RT * (adp - atp)
        dg_nadtrdh == dG0_nadtrdh + RT * (nad + nadph - nadh - nadp)
        dg_lp      == dG0_lp + RT * (lac_ext - lac)

        # flux bounds due to kinetics and thermodynamics
        v_pts == pts * kcat_pts * thermo_factor(dg_pts/RT)
        v_ed == ed * kcat_ed * thermo_factor(dg_ed/RT)
        v_pyk == pyk * kcat_pyk * thermo_factor(dg_pyk/RT)
        v_ldh == ldh * kcat_ldh * thermo_factor(dg_ldh/RT)
        v_ppc == ppc * kcat_ppc * thermo_factor(dg_ppc/RT)
        v_akgsyn == akgsyn * kcat_akgsyn * thermo_factor(dg_akgsyn/RT)
        v_gdhm == gdhm * kcat_gdhm * thermo_factor(dg_gdhm/RT)
        v_burn == burn * kcat_burn * thermo_factor(dg_burn/RT)
        v_nadtrdh == nadtrdh * kcat_nadtrdh * thermo_factor(dg_nadtrdh/RT)
        v_lp == lp * kcat_lp * thermo_factor(dg_lp/RT)

        # mass balance constraints
        v_pts - v_ed                                 ==  mu # g6p
        v_ed - v_pyk - v_ppc                         ==  0  # pep
        v_ed + v_pyk - v_ldh - v_akgsyn              ==  0  # pyr
        v_ldh - v_lp                                 ==  0  # lac
        v_ppc - v_akgsyn                             ==  0  # oac
        v_akgsyn - v_gdhm                            ==  0  # akg
        v_gdhm                                       ==  0  # glu
        v_ed + v_pyk - v_pts - v_burn                ==  mu # atp
        v_pts + v_burn - v_ed - v_pyk                == -mu # adp
        v_ldh + v_gdhm - v_ed - v_nadtrdh            ==  0  # nad
        v_ed + v_nadtrdh - v_ldh - v_gdhm            ==  0  # nadh
        v_gdhm + v_nadtrdh - v_ed                   ==  0  # nadp
        v_ed + v_nadtrdh                             ==  0  # nadph

        # density constraint(s) (can add more, e.g. membrane)
        pts + ed + pyk + ldh + ppc + akgsyn + gdhm + burn + nadtrdh + lp  <= total_proteome_mass_fraction

        # minimum maintenance
        min_burn_flux <= v_burn

        # measured ratio constraints (not strictly necessary)
        # atp   == log(10)  + adp  # atp/adp = 10, but remember that the concentration variables are logged
        # nadh  == log(0.2) + nad
        # nadph == log(0.2) + nadp
    end

    @objective(glu_ed, Max, mu)
    optimize!(glu_ed)
    objective_value(glu_ed)

    results_ed = Dict(
        "mu"        => value(mu),
        "pts"       => value(pts),
        "ed"       => value(ed),
        "pyk"       => value(pyk),
        "ldh"       => value(ldh),
        "ppc"       => value(ppc),
        "akgsyn"    => value(akgsyn),
        "gdhm"      => value(gdhm),
        "burn"      => value(burn),
        "nadtrdh"   => value(nadtrdh),
        "lp"        => value(lp),
        "glc_ext"   => value(glc_ext),
        "g6p"       => value(g6p),
        "pyr"       => value(pyr),
        "pep"       => value(pep),
        "oac"       => value(oac),
        "akg"       => value(akg),
        "lac"       => value(lac),
        "glu"       => value(glu),
        "atp"       => value(atp),
        "adp"       => value(adp),
        "nad"       => value(nad),
        "nadh"      => value(nadh),
        "nadp"      => value(nadp),
        "nadph"     => value(nadph),
        "v_pts"     => value(v_pts),
        "v_ed"      => value(v_ed),
        "v_pyk"     => value(v_pyk),
        "v_ldh"     => value(v_ldh),
        "v_ppc"     => value(v_ppc),
        "v_akgsyn"  => value(v_akgsyn),
        "v_gdhm"    => value(v_gdhm),
        "v_burn"    => value(v_burn),
        "v_nadtrdh" => value(v_nadtrdh),
        "v_lp"      => value(v_lp))
#     return results_ed
# end

# end #module