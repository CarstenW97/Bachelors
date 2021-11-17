for prot in prots
    println(prot, " ", value(m[prot]))
end

for dg in dgs
    println(dg, " ", value(m[dg]))
end

for met in mets
    println(met, " ", exp(value(m[met])) * 1000)
end

for flux in fluxes
    println(flux, " ", value(m[flux]))
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
        return "GLUDy"
    elseif s == :v_act
        return "ACt2r"
    elseif s == :v_ackr
        return "ACKr"
    elseif s == :v_ptar
        return "PTAr"
    elseif s == :v_acald
        return "ACALD"
    elseif s == :v_alcd
        return "ALCD2x"
    elseif s == :v_etoht
        return "ETOHt2r"
    elseif s == :v_thd
        return "THD2"
    elseif s == :v_pfl
        return "PFL"
    elseif s == :v_fort
        return "FORt2"
    else
        return uppercase(string(s)[3:end])
    end
end
