module Mos_Model
using JuMP
using KNITRO


function mos_model(glc_ext)
    model = Model(
        optimizer_with_attributes(KNITRO.Optimizer,
        "ms_enable" => 1,
        "honorbnds" => 1,
        "ms_maxsolves" => 10), # might need more resolves when model gets more complicated
    )

    # make functions
    ϵ = 1e-6
    fmax(x) = (x + sqrt(ϵ + x^2))/2 # max(0, x) approximation (this creates a normal Julia function)
    register(model, :fmax, 1, fmax; autodiff = true) # JuMP emits a warning telling you to do this
    fmin(x) = -fmax(-x) # min(0, x) approximation
    register(model, :fmin, 1, fmin; autodiff = true) # JuMP emits a warning telling you to do this

    # bounds
    lb_enz = 0.0 # [g enz / g DW cell]
    ub_enz = 1.0 # [g enz / g DW cell]
    lb_met = log(1e-6) # [log(1 uM)]
    ub_met = log(10e-3) # [log(10 mM)]
    lb_dg = -100 # [kJ/mol]
    ub_dg = 100 # [kJ/mol]
    lb_v = -50 # [mmol/gDW/h]
    ub_v = 50 # [mmol/gDW/h]

    @variables model begin
        # biomass [1/h] (g6p + ATP -> ADP)
        0 <= mu <= 10

        # Enzyme concentrations [g enz / g DW]
        lb_enz <= pts <= ub_enz # enzyme for: glc external + ATP -> ADP + g6p
        lb_enz <= emp <= ub_enz # enzyme for: g6p + 2NAD + 2ADP -> 2NADH + 2ATP + 2pep
        lb_enz <= pyk <= ub_enz # enzyme for: pep + ADP -> ATP + pyr
        lb_enz <= ldh <= ub_enz # enzyme for: pyr + NADH -> NAD + lac
        lb_enz <= burn <= ub_enz # enzyme for: ATP -> ADP
        lb_enz <= lp <= ub_enz # enzyme for: lac internal -> lac external (lactate permease)

        # Concentrations of internal metabolites [log(M)]
        lb_met <= g6p <= ub_met # glucose-6-phosphate
        lb_met <= pyr <= ub_met # pyruvate
        lb_met <= pep <= ub_met # phosphoenolpyrovate
        lb_met <= lac <= ub_met # lactate
        lb_met <= atp <= ub_met # ATP
        lb_met <= adp <= ub_met # ADP
        lb_met <= nad <= ub_met # NAD
        lb_met <= nadh <= ub_met # NADH

        lb_met <= glc_ext <= ub_met

        # Fluxes [mmol/gDW/h]
        lb_v <= v_pts <= ub_v
        lb_v <= v_emp <= ub_v
        lb_v <= v_pyk <= ub_v
        lb_v <= v_ldh <= ub_v
        lb_v <= v_burn <= ub_v
        lb_v <= v_lp <= ub_v

        # Thermodynamic variables
        lb_dg <= dg_pts <= ub_dg
        lb_dg <= dg_emp <= ub_dg
        lb_dg <= dg_pyk <= ub_dg
        lb_dg <= dg_ldh <= ub_dg
        lb_dg <= dg_burn <= ub_dg
        lb_dg <= dg_lp <= ub_dg
    end

    @NLparameters model begin
        # media conditions
        #glc_ext == log(50e-3) # [log(50 mM)]
        lac_ext == log(10e-6) # [log(10 uM)]

        # capacity constraint
        total_proteome_mass_fraction == 0.26 # [g enz/gDW]

        # maintenance requirement
        min_burn_flux == 1 # [mmol/gDW/h]

        # thermodynamic properties
        RT == 8.3145e-3 * 298.15 # gas constant [kJ/K/mol] * temperature [K]

        # enzyme specific activities (related to kcats) [mmol/g enz/h]
        kcat_pts == 213.75
        kcat_emp == 126
        kcat_pyk == 126
        kcat_ldh == 31
        kcat_burn == 10
        kcat_lp == 100

        # Standard Gibbs energy of reaction [kJ/mol]
        dg0_pts == -16.7
        dg0_emp == -4.71
        dg0_pyk == -31.7
        dg0_ldh == -23.7
        dg0_burn == -57
        dg0_lp == -10.3
    end

    @NLconstraints model begin
        # ΔG at concentrations
        dg_pts == dg0_pts + RT * (adp + g6p - atp - glc_ext)
        dg_emp == dg0_emp + RT * (2 * nadh + 2 * atp + 2 * pep - 2 * adp - 2 * nad - g6p)
        dg_pyk == dg0_pyk + RT * (pyr + atp - adp - pep)
        dg_ldh == dg0_ldh + RT * (lac + nad - nadh - pyr)
        dg_burn == dg0_burn + RT * (adp - atp)
        dg_lp == dg0_lp + RT * (lac_ext - lac)

        # flux bounds due to kinetics and thermodynamics
        pts * kcat_pts * fmin(exp(-dg_pts/RT) - 1) <= v_pts
        emp * kcat_emp * fmin(exp(-dg_emp/RT) - 1) <= v_emp
        pyk * kcat_pyk * fmin(exp(-dg_pyk/RT) - 1) <= v_pyk
        ldh * kcat_ldh * fmin(exp(-dg_ldh/RT) - 1) <= v_ldh
        burn * kcat_burn * fmin(exp(-dg_burn/RT) - 1) <= v_burn
        lp * kcat_lp * fmin(exp(-dg_lp/RT) - 1) <= v_lp

        v_pts <= pts * kcat_pts * fmax(1 - exp(dg_pts/RT))
        v_emp <= emp * kcat_emp * fmax(1 - exp(dg_emp/RT))
        v_pyk <= pyk * kcat_pyk * fmax(1 - exp(dg_pyk/RT))
        v_ldh <= ldh * kcat_ldh * fmax(1 - exp(dg_ldh/RT))
        v_burn <= burn * kcat_burn * fmax(1 - exp(dg_burn/RT))
        v_lp <= lp * kcat_lp * fmax(1 - exp(dg_lp/RT))

        # mass balance constraints
        v_pts - v_emp == mu # g6p
        2 * v_emp - v_pyk == 0 # pep
        v_pyk - v_ldh == 0 # pyr
        v_ldh - v_lp == 0 # lac
        2 * v_emp + v_pyk - v_pts - v_burn == mu # atp
        v_pts + v_burn - 2 * v_emp - v_pyk == -mu # adp
        v_ldh - 2 * v_emp == 0 # nad
        2 * v_emp  - v_ldh == 0 # nadh

        # density constraint(s) (can add more, e.g. membrane)
        pts + emp + pyk + ldh + lp + burn <= total_proteome_mass_fraction

        # minimum maintenance
        min_burn_flux <= v_burn

        # measured ratio constraints (not strictly necessary)
        atp == log(10) + adp # atp/adp = 10, but remember that the concentration variables are logged
        nadh == log(0.2) + nad
    end

    @objective(model, Max, mu)
    optimize!(model)
    objective_value(model)

    results_mo = Dict(
        "mu"        => value(mu),
        "glc_ext"   => value(glc_ext),
        "g6p"       => value(g6p),
        "pyr"       => value(pyr),
        "pep"       => value(pep),
        "lac"       => value(lac), 
        "atp"       => value(atp),
        "adp"       => value(adp),
        "nad"       => value(nad),
        "nadh"      => value(nadh),       
        "v_pts"     => value(v_pts),
        "v_emp"     => value(v_emp),
        "v_pyk"     => value(v_pyk),
        "v_ldh"     => value(v_ldh),
        "v_burn"    => value(v_burn),
        "v_lp"      => value(v_lp))
    return results_mo
end     

end #module