using Oceananigans
using Oceananigans.Units
using ClimaOcean
using SeawaterPolynomials.TEOS10
using FjordsSim

const FT = Oceananigans.defaults.FloatType

arch = GPU()
grid = ImmersedBoundaryGrid(
    joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_105to232_bathymetry_v3.nc"),
    arch,
    (7, 7, 7),
)
buoyancy = SeawaterBuoyancy(FT, equation_of_state=TEOS10EquationOfState(FT))
closure = (
    TKEDissipationVerticalDiffusivity(minimum_tke=7e-6),
    Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν=15, κ=10),
)
tracer_advection = (T=WENO(), S=WENO(), e=nothing, ϵ=nothing)
momentum_advection = WENOVectorInvariant(FT)
tracers = (:T, :S, :e, :ϵ)
initial_conditions = (T=5.0, S=33.0)
free_surface = SplitExplicitFreeSurface(grid, cfl=0.7)
coriolis = HydrostaticSphericalCoriolis(FT)
forcing = forcing_from_file(;
    grid=grid,
    filepath=joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_105to232_forcing_v2.nc"),
    tracers=tracers,
)
tbbc = top_bottom_boundary_conditions(;
        grid=grid,
        bottom_drag_coefficient=0.003,
    )
boundary_conditions = map(x -> FieldBoundaryConditions(;x...), tbbc)
# sobc = south_open_boundary_conditions()
# boundary_conditions = map(x -> FieldBoundaryConditions(;x...), recursive_merge(tbbc, sobc))
atmosphere = JRA55PrescribedAtmosphere(arch, FT;
    latitude=(58.98, 59.94),
    longitude=(10.18, 11.03),
    dir=joinpath(homedir(), "FjordsSim_data", "JRA55"),
)
downwelling_radiation = Radiation(arch, FT;
    ocean_emissivity=0.96,
    ocean_albedo=0.1
)
sea_ice = FreezingLimitedOceanTemperature()
biogeochemistry = nothing
results_dir = joinpath(homedir(), "FjordsSim_results", "oslofjord")
stop_time = 365days

simulation = coupled_hydrostatic_simulation(
    grid,
    buoyancy,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    initial_conditions,
    free_surface,
    coriolis,
    forcing,
    boundary_conditions,
    atmosphere,
    downwelling_radiation,
    sea_ice,
    biogeochemistry;
    results_dir,
    stop_time,
)

simulation.callbacks[:progress] = Callback(progress, TimeInterval(1hour))
ocean_sim = simulation.model.ocean
ocean_model = ocean_sim.model

prefix = joinpath(results_dir, "snapshots_ocean")
ocean_sim.output_writers[:ocean] = NetCDFWriter(
    ocean_model,
    (
        T=ocean_model.tracers.T,
        S=ocean_model.tracers.S,
        u=ocean_model.velocities.u,
        v=ocean_model.velocities.v,
    );
    filename = "$prefix",
    schedule = TimeInterval(1hour),
    overwrite_existing = true,
)

conjure_time_step_wizard!(simulation; cfl = 0.1, max_Δt = 3minutes, max_change = 1.01)
run!(simulation)
