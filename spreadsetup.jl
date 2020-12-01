basedir = @__DIR__

#' # Setup
#' 
#' First, load the required packages. Dates is a core julia package that
#' give us date/time handling, GeoData simplifies the loading of geospatial
#' raster files. It requires also loading NCDatasets.jl and ArchGDAL.jl to load
#' NetCDF and GeoTiff files, respectively.

using DimensionalData, GeoData, ArchGDAL, Dispersal
using Statistics, Dates, Plots, Unitful, NCDatasets, Setfield
using DimensionalData: setdims, rebuild, Between
using ColorSchemes, Colors
using CSV

#' ### Define simulation settings
#' 
#' We use DateTime unitsof months for the timestep:

timestep = Week(1)
lonmin = 113.3402
lonmax = 153.9523
latmin = -43.62234
latmin = -43.56234 # Fix for geodata tiff/netcdf load save errors
latmax = -10.65125
aus = Lon(Between(lonmin, lonmax)), Lat(Between(latmin, latmax))
usa = Lon(Between(-125.0, -50)), Lat(Between(4.5, 50))
tropics = Lon(Between(-6.0, 147)), Lat(Between(-16, 34.0))
incursionpoints = (
    Melbourne=(-37.805896, 144.959527),
    Mildura=(-34.219504, 142.130864),
    Coffs_Harbour=(-30.287245, 153.092991),
    Sydney=(-33.839943, 151.006101),
    Adelaide=(-34.901608, 138.601547),
    Port_Augusta=(-32.466201, 137.813850),
    Devonport=(-41.180545, 146.314887),
    Hobart=(-42.881742, 147.323879),
    Brisbane=(-27.436190, 152.990588),
    Cairns=(-16.937281, 145.747709),
    Seisia=(-10.928891, 142.411627),
    Perth=(-31.9505, 115.8605),
    Geraldton=(-28.778138, 114.615632),
    St_Hyacinth_QC = (45.627716, -72.956093),
    Tifton_GA = (31.445379, -83.513415),
    Gainsville_FL = (29.661792, -82.326272),
    Homestead_FL=(25.470084, -80.464693),
    Texas=(26.243250, -97.911040)
)

incursionpointstatekeys = (
    Melbourne=:Vic,
    Mildura=:Vic,
    Coffs_Harbour=:NSW,
    Sydney=:NSW,
    Adelaide=:SA,
    Port_Augusta=:SA,
    Devonport=:Tas,
    Hobart=:Tas,
    Brisbane=:QLD,
    Cairns=:QLD,
    Perth=:WA,
    Geraldton=:WA,
    St_Hyacinth_QC=:Canada,
    Tifton_GA=:USA,
    Gainsville_FL=:USA,
    Homestead_FL=:USA,
    Texas=:USA
)

#' ## Define a RuleSet
#' 
#' This will involve combining multiple dispersal componenents into a single
#' `RuleSet` object: population growth, local dispersal, Allee effects, and human
#' dispersal.
#' 
#' ### Climate driven population growth
#' 
#' Load the growthrates .16layer from netcdf. Make sure the mode.span is correct.

growthratespath = string(joinpath(basedir, "output/growthrates"))
gr_slices = map(readdir(growthratespath; join=true)) do path
    GDALarray(path; mappedcrs=EPSG(4326))[Band(1)] |>
        a->permutedims(a, (Lat, Lon))
end
growthtimespan = DateTime(2017,1):Month(1):DateTime(2017,12)
growthrates = cat(gr_slices...; dims=Ti(growthtimespan))
plot(growthrates[Ti(3), tropics...])

growthratesaus = growthrates[aus...]
growthratestropics = growthrates[tropics...]
growthratesusa = growthrates[usa...]

plot(growthratesusa[Ti(1:3:12)]; legend=:none, clims=(0.0, 0.15))
plot(mean(growthratesusa; dims=Ti); clims=(0, 0.15))

#' ### Define masking layers
#' 
#' The boolean mask lets the simulation know which cells should be ignored.
#' The missing mask can be used to mask maps before plotting.

# tropics
boolmasktropics = GeoData.boolmask(growthratestropics[Ti(1)])
missingmasktropics = GeoData.missingmask(growthratestropics[Ti(1)])
growthrates_zeromissingtropics = replace_missing(growthratestropics, 0.0)
growthmasktropics = rebuild(boolmasktropics, mean(growthrates_zeromissingtropics; dims=Ti)[Ti(1)] .* boolmasktropics .> 0)

# usa
boolmaskusa = GeoData.boolmask(growthratesusa[Ti(1)])
missingmaskusa = GeoData.missingmask(growthratesusa[Ti(1)])
growthrates_zeromissingusa = replace_missing(growthratesusa, 0.0)
growthmaskusa = rebuild(boolmaskusa, mean(growthrates_zeromissingusa; dims=Ti)[Ti(1)] .* boolmaskusa .> 0)

#aus
missingval(growthratesaus)
boolmaskaus = GeoData.boolmask(growthratesaus[Ti(1)])
missingmaskaus = GeoData.missingmask(growthratesaus[Ti(1)])
growthrates_zeromissingaus = replace_missing(growthratesaus, 0.0)
growthmaskaus = rebuild(boolmaskaus, mean(growthrates_zeromissingaus; dims=Ti)[Ti(1)] .* boolmaskaus .> 0)

growthmaskaus |> plot

# Check out growthrates and cropvalue arrays match
# Broken...
# @assert all(dims(states, Lat).val .≈ dims(growthrates, Lat).val)
# @assert all(dims(states, Lon).val .≈ dims(growthrates, Lon).val)

#' ## Define Rules
#' 
#' Create a `ExactLogisticGrowthMap` rule from the layer, here we use
#' unitful units for the layers' time dimension:
carrycap = 1e8
growth = LogisticGrowthMap{:population,:population}(
    auxkey=Val(:growthrates),
    carrycap=carrycap,
    timestep=Day(1),
);

#' ### Local dispersal
#' 
#' Local dispersal simulates natural dispersal of populations, according
#' to ability to fly or use other mechanisms.

λ = 0.004
radius = 2
@time hood = DispersalKernel{radius}(
    formulation=ExponentialKernel(λ),
    distancemethod=AreaToArea(30),
)
localdisp = InwardsPopulationDispersal{:population,:population}(hood)
log.(hood.kernel) |> heatmap

#' ### Allee effects
#' 
#' Allee effects specify minimum population required to sustain growth
#' within a cell. Populations below the `minfounders` threshold will be removed.
allee = AlleeExtinction{:population,:population}(minfounders=10000.0);

# TODO update this and test
# This is modifying the package version
import DynamicGrids: gridsize, inbounds
@inline Dispersal.applyrule!(data, rule::JumpDispersal{R,W}, state, index
) where {R,W} = begin
    # Ignore empty cells
    state > zero(state) || return state
    # Random dispersal events
    # rand() < rule.prob_threshold || return state
    # Randomly select spotting distance
    intspot = round(Int, rule.spotrange)
    rnge = -intspot:intspot
    jump = (rand(rnge), rand(rnge))
    jumpdest, is_inbounds = inbounds(jump .+ index, gridsize(data), RemoveOverflow())
    # Update spotted cell if it's on the grid
    if is_inbounds
      @inbounds data[W][jumpdest...] += state * rule.prob_threshold
    end
    state * (1 - rule.prob_threshold)
end

#' # Define Jump rule
jump = JumpDispersal{:population,:population}(
    spotrange=60.0,
    prob_threshold=0.05
)

#' ### Define a combined ruleset
dispersalruleset = Ruleset(
    Chain(localdisp, allee, growth), jump;
    timestep=timestep
);

#' ## Output
#' Define some color processors to use in live simuulations.
zerocolor = ARGB32(0.7)
maskcolor = ARGB32(0.2)
textconfig = if Sys.iswindows()
    TextConfig(font="arial", bcolor=maskcolor, fcolor=zerocolor)
else
    textconfig = TextConfig(font="cantarell", bcolor=maskcolor, fcolor=zerocolor)
end
oranges = ColorProcessor(ColorSchemes.Oranges_3, zerocolor, maskcolor, textconfig)
jet = ColorProcessor(ColorSchemes.jet, zerocolor, maskcolor, textconfig)
viridis = ColorProcessor(ColorSchemes.viridis, zerocolor, maskcolor, textconfig)
inferno = ColorProcessor(ColorSchemes.inferno, zerocolor, maskcolor, textconfig)
magma = ColorProcessor(ColorSchemes.magma, zerocolor, maskcolor, textconfig)
blues = ColorProcessor(ColorSchemes.Blues_3, zerocolor, maskcolor, textconfig)
algae = ColorProcessor(ColorSchemes.algae, zerocolor, maskcolor, textconfig)
cyclic = ColorProcessor(ColorSchemes.cyclic_grey_15_85_c0_n256, zerocolor, maskcolor, textconfig)
rainbow1 = ColorProcessor(ColorSchemes.rainbow1, zerocolor, maskcolor, textconfig)
wistia = ColorProcessor(ColorSchemes.Wistia, zerocolor, maskcolor, textconfig)
autumn = ColorProcessor(ColorSchemes.Wistia, zerocolor, maskcolor, textconfig)
