module GlnModel

using JuMP
using KNITRO

function gln_model(;
    glc_ext = 0.05,
    lac_ext = 1e-4,
    nh4_ext = 0.01,
    ac_ext = 1e-6,
    etoh_ext = 1e-6,
    proteome = 0.5,
    num_ms = 10,
    atp_adp_ratio = 10,
    nadh_nad_ratio = 0.2,
)

    gln_model = Model(
        optimizer_with_attributes(
            KNITRO.Optimizer,
            "ms_enable" => 1,
            "ms_maxsolves" => num_ms,
        )
    )

    # JuMP.set_silent(gln_model)

    # Thermodynamic function
    thermo_factor(x) = tanh(-0.7 * x) # x = ΔrG/RT
    register(gln_model, :thermo_factor, 1, thermo_factor; autodiff = true)

    # Bounds
    LB_enz = 0.0 # [g enz / g DW cell]
    UB_enz = proteome # [g enz / g DW cell]
    LB_met = log(1e-9) # [log(1 nM)]
    UB_met = log(100e-3) # [log(100 mM)]
    LB_v = -1 # [mmol/gDW/h]
    UB_v = 100 # [mmol/gDW/h]
    LB_dg = -80 # [kJ/mol]
    UB_dg = 1 # [kJ/mol]

    @variables gln_model begin
        0 <= mu <= UB_v # gln production rate

        # atp burn (excess)
        0 <= v_burn <= UB_v # mmol/gDW/h

        # Enzymes concentrations [g enz / g DW]
        LB_enz <= pts <= UB_enz
        LB_enz <= pgi <= UB_enz
        LB_enz <= pfk <= UB_enz
        LB_enz <= fba <= UB_enz
        LB_enz <= tpi <= UB_enz
        LB_enz <= gapd <= UB_enz
        LB_enz <= pgk <= UB_enz
        LB_enz <= pgm <= UB_enz
        LB_enz <= eno <= UB_enz
        LB_enz <= pyk <= UB_enz
        LB_enz <= ppc <= UB_enz
        LB_enz <= ldh <= UB_enz
        LB_enz <= pdh <= UB_enz
        LB_enz <= lact <= UB_enz
        LB_enz <= cs <= UB_enz
        LB_enz <= aconta <= UB_enz
        LB_enz <= acontb <= UB_enz
        LB_enz <= icdh <= UB_enz
        LB_enz <= gludy <= UB_enz
        LB_enz <= glns <= UB_enz
        LB_enz <= nh4t <= UB_enz
        LB_enz <= ptar <= UB_enz
        LB_enz <= ackr <= UB_enz
        LB_enz <= act <= UB_enz
        LB_enz <= acald <= UB_enz
        LB_enz <= alcd <= UB_enz
        LB_enz <= etoht <= UB_enz
        LB_enz <= thd <= UB_enz

        # Metabolite concentrations [log(M)]
        LB_met <= g6p <= UB_met
        LB_met <= f6p <= UB_met
        LB_met <= fdp <= UB_met
        LB_met <= dhap <= UB_met
        LB_met <= g3p <= UB_met
        LB_met <= dpg13 <= UB_met
        LB_met <= pg3 <= UB_met
        LB_met <= pg2 <= UB_met
        LB_met <= pep <= UB_met
        LB_met <= pyr <= UB_met
        LB_met <= lac <= UB_met
        LB_met <= accoa <=  UB_met
        LB_met <= oaa <= UB_met
        LB_met <= cit <= UB_met
        LB_met <= acon <= UB_met
        LB_met <= icit <= UB_met
        LB_met <= akg <= UB_met
        LB_met <= nh4 <= UB_met
        LB_met <= glu <= UB_met
        LB_met <= gln <= UB_met
        LB_met <= adp <= UB_met
        LB_met <= atp <= UB_met
        LB_met <= nad <= UB_met
        LB_met <= nadh <= UB_met
        LB_met <= acetald <= UB_met
        LB_met <= etoh <= UB_met
        LB_met <= actp <= UB_met
        LB_met <= ac <= UB_met

        # Fluxes [mmol/gDW/h]
        LB_v <= v_pts <= UB_v
        LB_v <= v_nh4t <= UB_v
        LB_v <= v_lact <= UB_v
        LB_v <= v_act <= UB_v
        LB_v <= v_etoht <= UB_v

        LB_v <= v_pgi <= UB_v
        LB_v <= v_pfk <= UB_v
        LB_v <= v_fba <= UB_v
        LB_v <= v_tpi <= UB_v
        LB_v <= v_gapd <= UB_v
        LB_v <= v_pgk <= UB_v
        LB_v <= v_pgm <= UB_v
        LB_v <= v_pyk <= UB_v
        LB_v <= v_eno <= UB_v
        LB_v <= v_ppc <= UB_v
        LB_v <= v_ldh <= UB_v
        LB_v <= v_pdh <= UB_v
        LB_v <= v_cs <= UB_v
        LB_v <= v_aconta <= UB_v
        LB_v <= v_acontb <= UB_v
        LB_v <= v_icdh <= UB_v
        LB_v <= v_gludy <= UB_v
        LB_v <= v_glns <= UB_v
        LB_v <= v_ptar <= UB_v
        LB_v <= v_ackr <= UB_v
        LB_v <= v_acald <= UB_v
        LB_v <= v_alcd <= UB_v
        LB_v <= v_thd <= UB_v

        # Thermodynamic variables
        LB_dg <= dg_pts <= UB_dg
        LB_dg <= dg_nh4t <= UB_dg
        LB_dg <= dg_lact <= UB_dg
        LB_dg <= dg_etoht <= UB_dg
        LB_dg <= dg_act <= UB_dg

        LB_dg <= dg_pgi <= UB_dg
        LB_dg <= dg_pfk <= UB_dg
        LB_dg <= dg_fba <= UB_dg
        LB_dg <= dg_tpi <= UB_dg
        LB_dg <= dg_gapd <= UB_dg
        LB_dg <= dg_pgk <= UB_dg
        LB_dg <= dg_pgm <= UB_dg
        LB_dg <= dg_pyk <= UB_dg
        LB_dg <= dg_eno <= UB_dg
        LB_dg <= dg_ppc <= UB_dg
        LB_dg <= dg_ldh <= UB_dg
        LB_dg <= dg_pdh <= UB_dg
        LB_dg <= dg_cs <= UB_dg
        LB_dg <= dg_aconta <= UB_dg
        LB_dg <= dg_acontb <= UB_dg
        LB_dg <= dg_icdh <= UB_dg
        LB_dg <= dg_gludy <= UB_dg
        LB_dg <= dg_glns <= UB_dg
        LB_dg <= dg_ptar <= UB_dg
        LB_dg <= dg_ackr <= UB_dg
        LB_dg <= dg_acald <= UB_dg
        LB_dg <= dg_alcd <= UB_dg
        LB_dg <= dg_thd <= UB_dg

    end

    @NLparameters gln_model begin
        RT == 8.3145e-3 * 298.15 # [kJ/mol]

        # Enzyme rates [mmol/g/h]
        kcat_pts == 7240.319618601774
        kcat_pgi == 7602.809491918478
        kcat_pfk == 29779.773325683316
        kcat_fba == 12070.317665674209
        kcat_tpi == 25881.04927925183
        kcat_gapd == 7096.06140643298
        kcat_pgk == 10031.027221032862
        kcat_pgm == 31309.78286506871
        kcat_pyk == 8271.195139265732
        kcat_eno == 5415.070665517126
        kcat_ppc == 1671.6634868719905
        kcat_ldh == 1604.7523387985657
        kcat_pdh == 7349.716396809779
        kcat_lact == 3963.1636080052076
        kcat_cs == 2361.1129054755233
        kcat_aconta == 1998.844393316627
        kcat_acontb == 1672.316435538672
        kcat_icdh == 865.4413532355705
        kcat_gludy == 13709.063214013708
        kcat_glns == 1708.433754802496
        kcat_nh4t == 6173.782075500597
        kcat_ptar == 23851.29662502135
        kcat_ackr == 26834.756823345528
        kcat_act == 7090.909090909091
        kcat_acald == 11454.251270331857
        kcat_alcd == 1622.4755081413768
        kcat_etoht == 7090.909090909091
        kcat_thd == 4528.47715441333

        # Gibbs energy of reaction
        dg0_pts == -2.9
        dg0_nh4t == 0.0 # diffusion
        dg0_etoht == 4.0 # proton symport
        dg0_act == 4.0 # proton symport
        dg0_lact == 4.0 # proton symport

        dg0_pgi == 0.9
        dg0_pfk == -14.7
        dg0_fba == 25.5
        dg0_tpi == 6.1
        dg0_gapd == 2.3
        dg0_pgk == -17.8
        dg0_pgm == -2.85
        dg0_pyk == -26.8
        dg0_eno == -4.3
        dg0_ppc == -38.7
        dg0_ldh == -23.15
        dg0_pdh == -35.9
        dg0_cs == -37.9
        dg0_aconta == -12.04
        dg0_acontb == 1.9
        dg0_icdh == 3.25
        dg0_gludy == -37.7
        dg0_glns == -16.1
        dg0_acald == 21.7
        dg0_ptar == 14.9
        dg0_ackr == -16.2
        dg0_alcd == -20.9
        dg0_thd == -1.4

        # media conditions log[M]
        glc_e == log(glc_ext)
        lac_e == log(lac_ext)
        co2 == log(1e-4)
        nh4_e == log(nh4_ext)
        etoh_e == log(etoh_ext)
        ac_e == log(ac_ext)

        # intracellular conditions log[M]
        phos == log(1e-3)
        nadph == log(1.2e-4)
        nadp == log(2.1e-6)
        coa == log(1.4e-3)
    end

    @NLconstraints gln_model begin
        # ΔG at concentrations
        dg_pts == dg0_pts + RT * (g6p + pyr - pep - glc_e)
        dg_nh4t == dg0_nh4t + RT * (nh4 - nh4_e)
        dg_lact == dg0_lact + RT * (lac_e - lac)
        dg_etoht == dg0_etoht + RT * (etoh_e - etoh)
        dg_act == dg0_act + RT * (ac_e - ac)

        dg_pgi == dg0_pgi + RT * (f6p - g6p)
        dg_pfk == dg0_pfk + RT * (fdp + adp - atp - f6p)
        dg_fba == dg0_fba + RT * (g3p + dhap - fdp)
        dg_tpi == dg0_tpi + RT * (g3p - dhap)
        dg_gapd == dg0_gapd + RT * (dpg13 + nadh - nad - phos - g3p)
        dg_pgk == dg0_pgk + RT * (pg3 + atp - adp - dpg13)
        dg_pgm == dg0_pgm + RT * (pg2 - pg3)
        dg_eno == dg0_eno + RT * (pep - pg2)
        dg_pyk == dg0_pyk + RT * (pyr + atp - adp - pep)
        dg_ppc == dg0_ppc + RT * (oaa + phos - co2 - pep)
        dg_ldh == dg0_ldh + RT * (lac + nad - nadh - pyr)
        dg_pdh == dg0_pdh + RT * (accoa + nadh + co2 - nad - coa - pyr)
        dg_cs == dg0_cs + RT * (cit + coa - oaa)
        dg_aconta == dg0_aconta + RT * (acon - cit)
        dg_acontb == dg0_acontb + RT * (icit - acon)
        dg_icdh == dg0_icdh + RT * (akg + nadph + co2 - nadp - icit)
        dg_gludy == dg0_gludy + RT * (glu + nadp - nadph - nh4 - akg)
        dg_glns == dg0_glns + RT * (gln + phos + adp - atp - nh4 - glu)
        dg_acald == dg0_acald + RT * (acetald + nad - nadh - accoa)
        dg_alcd == dg0_alcd + RT * (etoh + nad - nadh - acetald)
        dg_ptar == dg0_ptar + RT * (actp + coa - phos - accoa)
        dg_ackr == dg0_ackr + RT * (ac + atp - adp - actp)
        dg_thd == dg0_thd + RT * (nad + nadph - nadh - nadp)

        # flux to kinetics and thermodynamics
        v_pts == kcat_pts * pts * thermo_factor(dg_pts / RT)
        v_lact == kcat_lact * lact * thermo_factor(dg_lact / RT)
        v_nh4t == kcat_nh4t * nh4t * thermo_factor(dg_nh4t / RT)
        v_etoht == kcat_etoht * etoht * thermo_factor(dg_etoht / RT)
        v_act == kcat_act * act * thermo_factor(dg_act / RT)
        v_pgi == kcat_pgi * pgi * thermo_factor(dg_pgi / RT)
        v_pfk == kcat_pfk * pfk * thermo_factor(dg_pfk / RT)
        v_fba == kcat_fba * fba * thermo_factor(dg_fba / RT)
        v_tpi == kcat_tpi * tpi * thermo_factor(dg_tpi / RT)
        v_gapd == kcat_gapd * gapd * thermo_factor(dg_gapd / RT)
        v_pgk == kcat_pgk * pgk * thermo_factor(dg_pgk / RT)
        v_pgm == kcat_pgm * pgm * thermo_factor(dg_pgm / RT)
        v_eno == kcat_eno * eno * thermo_factor(dg_eno / RT)
        v_pyk == kcat_pyk * pyk * thermo_factor(dg_pyk / RT)
        v_ppc == kcat_ppc * ppc * thermo_factor(dg_ppc / RT)
        v_ldh == kcat_ldh * ldh * thermo_factor(dg_ldh / RT)
        v_pdh == kcat_pdh * pdh * thermo_factor(dg_pdh / RT)
        v_cs == kcat_cs * cs * thermo_factor(dg_cs / RT)
        v_aconta == kcat_aconta * aconta * thermo_factor(dg_aconta / RT)
        v_acontb == kcat_acontb * acontb * thermo_factor(dg_acontb / RT)
        v_icdh == kcat_icdh * icdh * thermo_factor(dg_icdh / RT)
        v_gludy == kcat_gludy * gludy * thermo_factor(dg_gludy / RT)
        v_acald == kcat_acald * acald * thermo_factor(dg_acald / RT)
        v_alcd == kcat_alcd * alcd * thermo_factor(dg_alcd / RT)
        v_ptar == kcat_ptar * ptar * thermo_factor(dg_ptar / RT)
        v_ackr == kcat_ackr * ackr * thermo_factor(dg_ackr / RT)
        v_thd == kcat_thd * thd * thermo_factor(dg_thd / RT)

        # mass balance constraints
        v_pts - v_pgi == 0 # g6p
        v_pgi - v_pfk == 0 # f6p
        v_pfk - v_fba == 0 # fdp
        v_fba - v_tpi == 0 # dhap
        v_fba + v_tpi - v_gapd == 0 # g3p
        v_gapd - v_pgk == 0 # dpg13
        v_pgk - v_pgm == 0 # pg3
        v_pgm - v_eno == 0 # pg2
        v_eno - v_pyk - v_ppc == 0 # pep
        v_pyk - v_ldh - v_pdh == 0 # pyr
        v_ldh - v_lact == 0 # lac
        v_pdh - v_cs - v_acald - v_ptar == 0 # accoa
        v_ppc - v_cs == 0 # oaa
        v_cs - v_aconta == 0 # cit
        v_aconta - v_acontb == 0 # acon
        v_acontb - v_icdh == 0 # icit
        v_icdh - v_gludy == 0 # akg
        v_nh4t - v_gludy - v_glns == 0 # nh4
        v_gludy - v_glns == 0 # glu
        v_glns - mu == 0 # gln
        v_pgk + v_pyk - v_pfk - v_glns - v_burn + v_ackr == 0 # atp
        -v_pgk - v_pyk + v_pfk + v_glns + v_burn - v_ackr == 0 # adp
        v_ldh - v_pdh - v_gapd + v_acald + v_alcd + v_thd == 0 # nad
        -v_ldh + v_pdh + v_gapd - v_acald - v_alcd - v_thd == 0 # nadh
        v_ptar - v_ackr == 0 # actp
        v_ackr - v_act == 0 # ac
        v_acald - v_alcd == 0 # acetald
        v_alcd - v_etoht == 0 # etoh

        # density constraint(s) (can add more, e.g. membrane)
        pts + pgi + pfk + fba + tpi + gapd + pgk + pgm + eno + pyk +
        ppc + ldh + pdh + lact + cs + aconta + acontb + icdh + gludy + glns +
        nh4t + ptar + ackr + acald + alcd + etoht + act + thd <= proteome

        # set ATP/ADP ratio NB: don't set the concentration bounds so that the ratio is not possible
        atp == log(atp_adp_ratio) + adp
        nadh == log(nadh_nad_ratio) + nad
    end

    @objective(gln_model, Max, mu)
    optimize!(gln_model)

    return gln_model
end

end #module
