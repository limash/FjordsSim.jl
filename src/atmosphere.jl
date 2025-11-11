using ClimaOcean.DataWrangling.JRA55: compute_bounding_nodes, infer_longitudinal_topology

import ClimaOcean.DataWrangling.JRA55: compute_bounding_indices

# Fix ClimaOcean for the custom longitude and latitude
# this is called from set! and uses grid to find the locations,
# which are 1 index more than necessary
function compute_bounding_indices(longitude::Nothing, latitude::Nothing, grid, LX, LY, λc, φc)
    λbounds = compute_bounding_nodes(longitude, grid, LX, λnodes)
    φbounds = compute_bounding_nodes(latitude, grid, LY, φnodes)

    i₁, i₂ = compute_bounding_indices(λbounds, λc)
    j₁, j₂ = compute_bounding_indices(φbounds, φc)
    TX = infer_longitudinal_topology(λbounds)

    # to prevent taking larger than grid areas
    i₁ = (i₂ - i₁ >= grid.Nx) ? (i₂ - grid.Nx + 1) : i₁
    j₁ = (j₂ - j₁ >= grid.Ny) ? (j₂ - grid.Ny + 1) : j₁

    return i₁, i₂, j₁, j₂, TX
end
