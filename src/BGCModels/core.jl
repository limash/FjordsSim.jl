"""
OxyDep basic biogeochemical transformations between NUT, P, HET, DOM, POM, O2
"""
# Limiting equations and switches
@inline yy(value, consta) = consta^2 / (value^2 + consta^2)   #This is a squared Michaelis-Menten type of limiter
@inline F_ox(conc, threshold) = (0.5 + 0.5 * tanh(conc - threshold))
@inline F_subox(conc, threshold) = (0.5 - 0.5 * tanh(conc - threshold))

# P
@inline LimLight(PAR, Iopt) = PAR / Iopt * exp(1.0 - PAR / Iopt)  #!Dependence of P growth on Light (Steel)
@inline LimN(Knut, NUT, P) = yy(Knut, NUT / max(0.0001, P)) #!Dependence of P growth on NUT
@inline Q₁₀(T) = 1.88^(T / 10) # T in °C  # inital for NPZD
#@inline LimT(T) = max(0., 2^((T-10.0)/10.) - 2^((T-32.)/3.)) # ERSEM
# = q10^((T-t_upt_min)/10)-q10^((T-t_upt_max)/3):  q10=2. !Coefficient for uptake rate dependence on t
# t_upt_min=10. !Low  t limit for uptake rate dependence on t; t_upt_max=32 !High t limit for uptake rate dependence on t
@inline LimT(T) = exp(0.0663 * (T - 0.0)) #for Arctic (Moore et al.,2002; Jin et al.,2008) 
# = exp(temp_aug_rate*(T-t_0)):  t_0= 0. !reference temperature temp_aug_rate = 0.0663 !temperature augmentation rate
#@inline light_limitation(PAR, α, Max_uptake) = α * PAR / sqrt(Max_uptake ^ 2 + α ^ 2 * PAR ^ 2)

#@inline GrowthPhy(Max_uptake,PAR,α,T,Knut,NUT,P,Iopt) = Max_uptake*LimT(T)*LimN(Knut,NUT,P)*light_limitation(PAR,α,Max_uptake)*P*Iopt/Iopt
@inline GrowthPhy(Max_uptake, PAR, α, T, Knut, NUT, P, Iopt) =
    Max_uptake * LimT(T) * LimN(Knut, NUT, P) * LimLight(PAR, Iopt) * α / α
@inline RespPhy(r_phy_nut, P) = r_phy_nut * P
@inline MortPhy(r_phy_pom, P) = r_phy_pom * P
@inline ExcrPhy(r_phy_dom, P) = r_phy_dom * P

# HET
@inline GrazPhy(r_phy_het, Kphy, P, HET) =
    r_phy_het * yy(Kphy, max(0.0, P - 0.01) / max(0.0001, HET)) * HET
@inline GrazPOM(r_pom_het, Kpom, POM, HET) =
    r_pom_het * yy(Kpom, max(0.0, POM - 0.01) / max(0.0001, HET)) * HET
@inline RespHet(r_het_nut, HET) = r_het_nut * HET
@inline MortHet(r_het_pom, HET, O₂, O2_suboxic) =
    (r_het_pom + F_subox(O₂, O2_suboxic) * 0.01 * r_het_pom) * HET

# POM
@inline POM_decay_ox(r_pom_nut_oxy, POM) = r_pom_nut_oxy * POM
@inline POM_decay_denitr(r_pom_nut_nut, POM, O₂, O2_suboxic, NUT) =
    r_pom_nut_nut * POM * F_subox(O₂, O2_suboxic) * F_ox(NUT, 0.01)
#! depends on NUT (NO3+NO2) and DOM (NH4+Urea+"real"DON) ! depends on T ! stops at NUT<0.01 
@inline Autolys(r_pom_dom, POM) = r_pom_dom * POM

# DOM
@inline DOM_decay_ox(r_dom_nut_oxy, DOM) = r_dom_nut_oxy * DOM
@inline DOM_decay_denitr(r_dom_nut_nut, DOM, O₂, O2_suboxic, NUT) =
    r_dom_nut_nut * DOM * F_subox(O₂, O2_suboxic) * F_ox(NUT, 0.01)
#! depends on NUT (NO3+NO2) and DOM (NH4+Urea+"real"DON) ! depends on T ! stops at NUT<0.01 

# O₂

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

@inline function (bgc::OXYDEP)(::Val{:NUT}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, PAR)
    Max_uptake = bgc.Max_uptake
    Knut = bgc.Knut
    α = bgc.initial_photosynthetic_slope
    r_phy_nut = bgc.r_phy_nut
    r_het_nut = bgc.r_het_nut
    r_pom_nut_oxy = bgc.r_pom_nut_oxy
    r_dom_nut_oxy = bgc.r_dom_nut_oxy
    NtoN = bgc.NtoN
    r_pom_nut_nut = bgc.r_pom_nut_nut
    O2_suboxic = bgc.O2_suboxic
    r_dom_nut_nut = bgc.r_dom_nut_nut
    Iopt = bgc.Iopt

    #println(GrowthPhy(Max_uptake,PAR,α,T,Knut,NUT,P,Iopt))
    #wait_for_key("press any key to continue")

    return (
        RespPhy(r_phy_nut, P) +
        RespHet(r_het_nut, HET) +
        DOM_decay_ox(r_dom_nut_oxy, DOM) +
        POM_decay_ox(r_pom_nut_oxy, POM) - GrowthPhy(Max_uptake, PAR, α, T, Knut, NUT, P, Iopt) -
        NtoN * (
            POM_decay_denitr(r_pom_nut_nut, POM, O₂, O2_suboxic, NUT) +
            DOM_decay_denitr(r_dom_nut_nut, DOM, O₂, O2_suboxic, NUT)
        )
    )
    # Denitrification of POM and DOM leads to decrease of NUT (i.e. NOx)
end

@inline function (bgc::OXYDEP)(::Val{:P}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, PAR)
    Max_uptake = bgc.Max_uptake
    Knut = bgc.Knut
    α = bgc.initial_photosynthetic_slope
    r_phy_het = bgc.r_phy_het
    Kphy = bgc.Kphy
    r_phy_nut = bgc.r_phy_nut
    r_phy_pom = bgc.r_phy_pom
    r_phy_dom = bgc.r_phy_dom
    Iopt = bgc.Iopt

    return (
        GrowthPhy(Max_uptake, PAR, α, T, Knut, NUT, P, Iopt) -
        GrazPhy(r_phy_het, Kphy, P, HET) - RespPhy(r_phy_nut, P) - MortPhy(r_phy_pom, P) -
        ExcrPhy(r_phy_dom, P)
    )
end

@inline function (bgc::OXYDEP)(::Val{:HET}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, PAR)
    r_phy_het = bgc.r_phy_het
    Kphy = bgc.Kphy
    r_pom_het = bgc.r_pom_het
    Kpom = bgc.Kpom
    r_het_nut = bgc.r_het_nut
    r_het_pom = bgc.r_het_pom
    Uz = bgc.Uz
    O2_suboxic = bgc.O2_suboxic

    return (
        Uz * (GrazPhy(r_phy_het, Kphy, P, HET) + GrazPOM(r_pom_het, Kpom, POM, HET)) -
        MortHet(r_het_pom, HET, O₂, O2_suboxic) - RespHet(r_het_nut, HET)
    )
end

@inline function (bgc::OXYDEP)(::Val{:POM}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, PAR)
    r_phy_het = bgc.r_phy_het
    Kphy = bgc.Kphy
    r_pom_het = bgc.r_pom_het
    Kpom = bgc.Kpom
    Uz = bgc.Uz
    Hz = bgc.Hz
    r_phy_pom = bgc.r_phy_pom
    r_het_pom = bgc.r_het_pom
    r_pom_nut_oxy = bgc.r_pom_nut_oxy
    r_pom_dom = bgc.r_pom_dom
    r_pom_nut_nut = bgc.r_pom_nut_nut
    O2_suboxic = bgc.O2_suboxic

    return (
        (1.0 - Uz) *
        (1.0 - Hz) *
        (GrazPhy(r_phy_het, Kphy, P, HET) + GrazPOM(r_pom_het, Kpom, POM, HET)) +
        MortPhy(r_phy_pom, P) +
        MortHet(r_het_pom, HET, O₂, O2_suboxic) - POM_decay_ox(r_pom_nut_oxy, POM) -
        Autolys(r_pom_dom, POM) - GrazPOM(r_pom_het, Kpom, POM, HET) -
        POM_decay_denitr(r_pom_nut_nut, POM, O₂, O2_suboxic, NUT)
    )
end

@inline function (bgc::OXYDEP)(::Val{:DOM}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, PAR)
    r_phy_het = bgc.r_phy_het
    Kphy = bgc.Kphy
    r_pom_het = bgc.r_pom_het
    Kpom = bgc.Kpom
    Uz = bgc.Uz
    Hz = bgc.Hz
    r_phy_dom = bgc.r_phy_dom
    r_dom_nut_oxy = bgc.r_dom_nut_oxy
    r_pom_dom = bgc.r_pom_dom
    r_pom_nut_nut = bgc.r_pom_nut_nut
    O2_suboxic = bgc.O2_suboxic

    return (
        (1.0 - Uz) *
        Hz *
        (GrazPhy(r_phy_het, Kphy, P, HET) + GrazPOM(r_pom_het, Kpom, POM, HET)) +
        ExcrPhy(r_phy_dom, P) - DOM_decay_ox(r_dom_nut_oxy, DOM) +
        Autolys(r_pom_dom, POM) +
        POM_decay_denitr(r_pom_nut_nut, POM, O₂, O2_suboxic, NUT)
    )
    # Denitrification of "real DOM" into NH4 (DOM_decay_denitr) will not change state variable DOM
end

@inline function (bgc::OXYDEP)(::Val{:O₂}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, PAR)
    Max_uptake = bgc.Max_uptake
    Knut = bgc.Knut
    α = bgc.initial_photosynthetic_slope
    r_phy_nut = bgc.r_phy_nut
    r_het_nut = bgc.r_het_nut
    r_pom_nut_oxy = bgc.r_pom_nut_oxy
    r_dom_nut_oxy = bgc.r_dom_nut_oxy
    OtoN = bgc.OtoN
    O2_suboxic = bgc.O2_suboxic
    Iopt = bgc.Iopt

    return (
        -OtoN * (
            RespPhy(r_phy_nut, P) +
            RespHet(r_het_nut, HET) +
            DOM_decay_ox(r_dom_nut_oxy, DOM) +
            POM_decay_ox(r_pom_nut_oxy, POM) -
            GrowthPhy(Max_uptake, PAR, α, T, Knut, NUT, P, Iopt) # due to OM production and decay in normoxia
            +
            DOM_decay_ox(r_dom_nut_oxy, DOM) * (F_subox(O₂, O2_suboxic))
        )
    )
    # (POM_decay_denitr + DOM_decay_denitr) & !denitrification doesn't change oxygen
    # (DOM_decay_ox(r_dom_nut_oxy,DOM)*(F_subox) !additional consumption of O₂ due to oxidation of reduced froms of S,Mn,Fe etc.
    # in suboxic conditions (F_subox) equals consumption for NH4 oxidation (Yakushev et al, 2008)

end

############################################################

"""
@inline function (bgc::OXYDEP)(::Val{:Ch_free}, x, y, z, t, NUT, P, HET, POM, DOM, O₂, T, Ch_free, PAR)
    Max_uptake = bgc.Max_uptake
    Knut = bgc.Knut
    α = bgc.initial_photosynthetic_slope
    r_phy_het = bgc.r_phy_het
    Kphy = bgc.Kphy
    r_phy_nut = bgc.r_phy_nut
    r_phy_pom = bgc.r_phy_pom
    r_phy_dom = bgc.r_phy_dom
    Iopt = bgc.Iopt

    return (
        0.1 * GrowthPhy(Max_uptake, PAR, α, T, Knut, NUT, P, Iopt) -
        GrazPhy(r_phy_het, Kphy, P, HET) - RespPhy(r_phy_nut, P) - MortPhy(r_phy_pom, P) -
        ExcrPhy(r_phy_dom, P)
    )
end
"""
