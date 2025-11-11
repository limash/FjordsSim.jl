module FjordsSim 

export
    # oceananigans methods
    ImmersedBoundaryGrid,
    # forcings
    forcing_from_file,
    # boundary conditions
    top_bottom_boundary_conditions,
    # simulations
    coupled_hydrostatic_simulation,
    # utils
    recursive_merge, progress

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Utils
using ClimaOcean
using NCDatasets
using Adapt

import Oceananigans.Advection: cell_advection_timescale

# to allow time step adjusting in OceanSeaIceModel
cell_advection_timescale(model::OceanSeaIceModel) = cell_advection_timescale(model.ocean.model)

include("atmosphere.jl")
include("boundary_conditions.jl")
include("forcing.jl")
include("grid.jl")
include("turbulence.jl")
include("utils.jl")

function coupled_hydrostatic_simulation(
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
    results_dir=joinpath(homedir(), "FjordsSim_results"),
    stop_time=365days,
)
    isdir(results_dir) || mkpath(results_dir)

    println("Start compiling HydrostaticFreeSurfaceModel")
    ocean_model = HydrostaticFreeSurfaceModel(;
        grid,
        buoyancy,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        free_surface,
        coriolis,
        forcing,
        boundary_conditions,
        biogeochemistry,
    )
    println("Done compiling HydrostaticFreeSurfaceModel")
    set!(ocean_model; initial_conditions...)
    Δt = 1second
    ocean_sim = Simulation(ocean_model; Δt, stop_time)
    interfaces = ComponentInterfaces(atmosphere, ocean_sim, sea_ice; radiation=downwelling_radiation)
    coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; atmosphere, radiation=downwelling_radiation, interfaces)
    println("Initialized coupled model")
    coupled_simulation = Simulation(coupled_model; Δt, stop_time)
    return coupled_simulation
end  # function coupled_hydrostatic_simulation

end  # module
