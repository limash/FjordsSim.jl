using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Oceananigans.Utils: prettytimeunits, maybe_int

include("setup.jl")

folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots")
S = FieldTimeSeries("$filename.jld2", "S")

salt = interior(S, :, :, :, :)

sim_setup = setup_region_3d()
grid = sim_setup.grid_callable(sim_setup.grid_args...)
