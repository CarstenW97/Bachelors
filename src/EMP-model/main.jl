using JuMP
using JSON

include(joinpath("src", "EMP-model", "model.jl"))
import .GlnModel

m = GlnModel.gln_model(
    glc_ext = 0.005,
    lac_ext = 1e-2,
    nh4_ext = 0.01,
    ac_ext = 1e-16,
    etoh_ext = 1e-6,
    atp_adp_ratio = 10,
    nadh_nad_ratio = 35, # around 35
)

prots = [:pts, :pgi, :pfk, :fba, :tpi, :gapd, :pgk, :pgm, :eno, :pyk, :ppc, :ldh, :pdh, :lact, :cs, :aconta, :acontb, :icdh, :gludy, :glns, :nh4t, :acald, :alcd, :etoht, :ackr, :ptar, :act, :thd]
mets = [:g6p, :f6p, :fdp, :dhap, :g3p, :dpg13, :pg3, :pg2, :pep, :pyr, :lac, :accoa, :oaa, :cit, :acon, :icit, :akg, :nh4, :glu, :gln, :acetald, :etoh, :actp, :ac, :adp, :atp, :nadh, :nad]
fluxes = [:v_pts, :v_pgi, :v_pfk, :v_fba, :v_tpi, :v_gapd,:v_pgk, :v_pgm, :v_pyk, :v_ppc, :v_ldh, :v_pdh, :v_lact, :v_cs, :v_aconta, :v_acontb, :v_icdh, :v_gludy, :v_glns, :v_nh4t, :v_eno, :v_acald, :v_thd, :v_alcd, :v_etoht, :v_ackr, :v_ptar, :v_act]
dgs = [:dg_pts, :dg_pgi, :dg_pfk, :dg_fba, :dg_tpi, :dg_gapd,:dg_pgk, :dg_pgm, :dg_pyk, :dg_ppc, :dg_ldh, :dg_pdh, :dg_lact, :dg_cs, :dg_aconta, :dg_acontb, :dg_icdh, :dg_gludy, :dg_glns, :dg_nh4t, :dg_thd, :dg_eno, :dg_acald, :dg_alcd, :dg_etoht, :dg_ackr, :dg_ptar, :dg_act]

for prot in prots
    println(prot, " ", value(m[prot]))
end

for dg in dgs
    println(dg, " ", value(m[dg]))
end

for met in mets
    println(met, " ", exp(value(m[met]))*1000)
end

for flux in fluxes
    println(flux, " ", value(m[fluxes]))
end

function map_rids(s)
    if s == :v_ldh
        return "LDH_D"
    elseif s == :v_lact
        return "D_LACt2"
    elseif s == :v_nh4t
        return "NH4t"
    elseif s == :v_pts
        return "GLCpts"
    elseif s == :v_aconta
        return "ACONTa"
    elseif s == :v_acontb
        return "ACONTb"
    elseif s == :v_icdh
        return "ICDHyr"
    elseif s == :v_gludy
        return  "GLUDy"
    elseif s == :v_act
        return  "ACt2r"
    elseif s == :v_ackr
        return  "ACKr"
    elseif s == :v_ptar
        return  "PTAr"
    elseif s == :v_acald
        return  "ACALD"
    elseif s == :v_alcd
        return  "ALCD2x"
    elseif s == :v_etoht
        return  "ETOHt2r"
    elseif s == :v_thd
        return  "THD2"
    else
        return uppercase(string(s)[3:end])
    end
end

open(joinpath("src", "EMP-model", "fluxes.json"), "w") do io
    JSON.print(io, Dict(
        map_rids(k) => value(m[k]) for k in fluxes
    ))
end
