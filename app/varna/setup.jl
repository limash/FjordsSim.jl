using Oceananigans.Architectures
using Oceananigans.Units

args_grid = (
    arch = GPU(),
    Nx = 119,
    Ny = 42,
    Nz = 20,
    dx = 200,
    dy = 50,
    z_levels = -reverse([
        0.0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
    ]),
    z_middle = -reverse([
        0.5,
        1.5,
        2.5,
        3.5,
        4.5,
        5.5,
        6.5,
        7.5,
        8.5,
        9.5,
        10.5,
        11.5,
        12.5,
        13.5,
        14.5,
        15.5,
        16.5,
        17.5,
        18.5,
        19.5,
    ]),
    halo = (7, 7, 7),
    datadir = joinpath(homedir(), "data_Varna"),
    filename = "Varna_topo.jld2",
)

args_oxydep = (
    initial_photosynthetic_slope = 0.1953 / day, # 1/(W/m²)/s
    Iopt = 50.0,     # (W/m2)
    alphaI = 1.8,   # [d-1/(W/m2)]
    betaI = 5.2e-4, # [d-1/(W/m2)]
    gammaD = 0.71,  # (-)
    Max_uptake = 2.5 / day,  # 1/d 2.0 4 5
    Knut = 2.0,            # (nd)
    r_phy_nut = 0.10 / day, # 1/d
    r_phy_pom = 0.15 / day, # 1/d
    r_phy_dom = 0.17 / day, # 1/d
    r_phy_het = 2.0 / day,  # 1/d 0.4
    Kphy = 0.1,             # (nd) 0.7
    r_pom_het = 0.7 / day,  # 1/d 0.7
    Kpom = 2.0,     # (nd)
    Uz = 0.6,       # (nd)
    Hz = 0.5,       # (nd)
    r_het_nut = 0.15 / day,      # 1/d 0.05
    r_het_pom = 0.10 / day,      # 1/d 0.02
    r_pom_nut_oxy = 0.006 / day, # 1/d
    r_pom_dom = 0.01 / day,      # 1/d
    r_dom_nut_oxy = 0.050 / day,  # 1/d
    O2_suboxic = 30.0,    # mmol/m3
    r_pom_nut_nut = 0.010 / day, # 1/d
    r_dom_nut_nut = 0.003 / day, # 1/d
    OtoN = 8.625, # (nd)
    CtoN = 6.625, # (nd)
    NtoN = 5.3,   # (nd)
    NtoB = 0.016, # (nd)
    sinking_speeds = (PHY = 0.15 / day, HET = 0.4 / day, POM = 10.0 / day),
)
