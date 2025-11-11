# FjordsSim.jl

## Installation

1. Clone the git repository `git clone https://github.com/AquaBadgerForge/FjordsSim.jl.git`.
2. Move to the downloaded folder `cd FjordsSim.jl`.
3. Run Julia REPL and activate the FjordsSim environment `julia --project=.`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.

## Run an example Oslofjord simulation

1. Download the [grid, forcing, atmospheric forcing](https://www.dropbox.com/scl/fo/gc3yc155b5eohi7998wgh/AGN2Yt3HyQ0LlZGImpcca6o?rlkey=x6okc3uxe2avud6sbxgd00l14&st=093llyqp&dl=0).
2. In `FjordsSim.jl/app/oslofjord.jl` it is possible to specify the location of the input data files.
By default, the files should be in `$HOME/FjordsSim_data/oslofjord/` and `$HOME/FjordsSim_data/JRA55/`.
Also, it is possible to specify the result folder destination.
By default, the result will go to `$HOME/FjorsSim_results/oslofjord/`.
3. Run `julia --project app/oslofjord.jl`.
This will generate a netcdf results files.

![example_result](./artifacts/phytoplankton_multi.png)
