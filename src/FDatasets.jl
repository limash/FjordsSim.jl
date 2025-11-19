module FDatasets

export DSForcing

using NCDatasets
using Oceananigans
using ClimaOcean.DataWrangling: Metadatum, metadata_path

import Oceananigans: location
import ClimaOcean.DataWrangling:
    metadata_filename,
    default_download_directory,
    all_dates,
    first_date,
    last_date,
    dataset_variable_name,
    download_dataset,
    longitude_interfaces,
    latitude_interfaces,
    z_interfaces,
    reversed_vertical_axis,
    inpainted_metadata_path

struct DSForcing
    metadata_filename::String
    default_download_directory::String
    reversed_vertical_axis::Any
    longitude_interfaces::Any
    latitude_interfaces::Any
    size::Any
    all_dates::Any
    first_date::Any
    last_date::Any
    z_interfaces::Any
end

function DSForcing(metadata_filename, default_download_directory)
    reversed_vertical_axis = false
    filepath = joinpath(default_download_directory, metadata_filename)
    ds = NCDataset(filepath)
    longitude_interfaces = (ds["Nx"][1], ds["Nx"][end])
    latitude_interfaces = (ds["Ny"][1], ds["Ny"][end])
    size = (ds.dim["Nx"], ds.dim["Ny"], ds.dim["Nz"])
    all_dates = ds["time"][:]
    first_date = ds["time"][1]
    last_date = ds["time"][end]
    z_interfaces = ds["Nz_faces"][:]
    close(ds)

    return DSForcing(
        metadata_filename,
        default_download_directory,
        reversed_vertical_axis,
        longitude_interfaces,
        latitude_interfaces,
        size,
        all_dates,
        first_date,
        last_date,
        z_interfaces,
    )
end  # function

Forcing_dataset_variable_names = Dict(:temperature => "T", :salinity => "S", :u_velocity => "u", :v_velocity => "v")
Variable_location = Dict(
    :temperature => (Center, Center, Center),
    :salinity => (Center, Center, Center),
    :free_surface => (Center, Center, Nothing),
    :sea_ice_thickness => (Center, Center, Nothing),
    :sea_ice_concentration => (Center, Center, Nothing),
    :net_heat_flux => (Center, Center, Nothing),
    :u_velocity => (Face, Center, Center),
    :v_velocity => (Center, Face, Center),
    :sensible_heat_flux => (Center, Center, Nothing),
    :latent_heat_flux => (Center, Center, Nothing),
    :net_longwave => (Center, Center, Nothing),
    :downwelling_shortwave => (Center, Center, Nothing),
    :downwelling_longwave => (Center, Center, Nothing),
)

const MetadatumForcing = Metadatum{<:DSForcing,<:Any,<:Any}

metadata_filename(metadatum::MetadatumForcing) = metadatum.dataset.metadata_filename
default_download_directory(ds::DSForcing) = ds.default_download_directory
reversed_vertical_axis(ds::DSForcing) = ds.reversed_vertical_axis
longitude_interfaces(ds::DSForcing) = ds.longitude_interfaces
latitude_interfaces(ds::DSForcing) = ds.latitude_interfaces
Base.size(ds::DSForcing) = ds.size
Base.size(ds::DSForcing, variable) = size(ds)

all_dates(ds::DSForcing, args...) = ds.all_dates
first_date(ds::DSForcing, args...) = ds.first_date
last_date(ds::DSForcing, args...) = ds.last_date

z_interfaces(metadatum::MetadatumForcing) = metadatum.dataset.z_interfaces

dataset_variable_name(metadatum::MetadatumForcing) = Forcing_dataset_variable_names[metadatum.name]
location(metadatum::MetadatumForcing) = Variable_location[metadatum.name]

function download_dataset(metadatum::MetadatumForcing)
    filepath = metadata_path(metadatum)
    return filepath
end

function inpainted_metadata_filename(metadatum::MetadatumForcing)
    original_filename = metadata_filename(metadatum)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadatum::MetadatumForcing) = joinpath(metadatum.dir, inpainted_metadata_filename(metadatum))

end  # module
