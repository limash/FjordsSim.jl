# Copyright 2024 The FjordsSim Authors.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using Oceananigans.Units: second, seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Utils: TimeInterval, IterationInterval
using Oceananigans.Simulations: Callback, conjure_time_step_wizard!, run!
using Oceananigans.OutputWriters: NetCDFWriter
using Oceanostics
using FjordsSim: coupled_hydrostatic_simulation, progress
using Printf

include("setup.jl")

## Model Setup
sim_setup = setup_region_3d_45to82()

coupled_simulation = coupled_hydrostatic_simulation(sim_setup)

coupled_simulation.callbacks[:progress] = Callback(progress, TimeInterval(3hours))

ocean_sim = coupled_simulation.model.ocean
ocean_model = ocean_sim.model

net_ocean_fluxes = NamedTuple((
    u_atm_ocean_flux = coupled_simulation.model.interfaces.net_fluxes.ocean_surface.u,
    v_atm_ocean_flux = coupled_simulation.model.interfaces.net_fluxes.ocean_surface.v,
))
prefix = joinpath(sim_setup.results_dir, "snapshots_ocean")
ocean_sim.output_writers[:ocean] = NetCDFWriter(
    ocean_model,
    merge(
        ocean_model.tracers,
        ocean_model.velocities,
        coupled_simulation.model.interfaces.atmosphere_ocean_interface.fluxes,
        net_ocean_fluxes,
    );
    filename = "$prefix",
    schedule = TimeInterval(1hours),
    overwrite_existing = true,
)

## Spinning up the simulation
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.1.
ocean_sim.stop_time = 10days
coupled_simulation.stop_time = 10days

conjure_time_step_wizard!(coupled_simulation; cfl = 0.1, max_Δt = 1minute, max_change = 1.01)
run!(coupled_simulation)

## Running the simulation
# This time, we set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))
ocean_sim.stop_time = 355days
coupled_simulation.stop_time = 355days

conjure_time_step_wizard!(coupled_simulation; cfl = 0.25, max_Δt = 3minutes, max_change = 1.01)
run!(coupled_simulation)
