basedir = @__DIR__
year = 2016

#' 
#' # Growth rate calculation
#' 
#' In this example we will calculate the expected population growth rates of our species.
#' 
#' Use biological data from Early 2017
#' GBIF distribution data
#' Spread data from Early 2017
#' 
#' ## Load some required packages
#' 
#' These packages take care of loading and plotting data, and handling sci units and dates.

using GrowthMaps, Plots, Unitful, UnitfulRecipes, Dates, Setfield, Statistics
using GeoData, ArchGDAL, NCDatasets, HDF5, CSV
using Unitful: °C, K, cal, mol

#' 
#' ## Define model components
#' 
#' First we'll define the growth model using `SchoolfieldIntrinsicGrowth`, based on
#' Schoolfield (1981).
#' 
#' When defining model components, the first parameter is a `:symbol` for the
#' required raster layer in the source data.

p27 = Param(0.17; bounds=(0.0, 1.0))
# ΔH_A = Param(15588.0; units=cal/mol, bounds=(1e5, 2e5))
# ΔH_H = Param(73744.0; units=cal/mol, bounds=(0.5, 1e5))
Thalf_H = Param(307.6; units=K, bounds=(290.0, 320.0))
# ΔH_L = Param(73744.0; units=cal/mol, bounds=(0.5e5, 1e5))
Thalf_L = Param(500.0; units=K, bounds=(100.0, 1000.0))
T_ref = K(27.0°C)

# Barfield 1978 (need to convert parameters using Schoolfield)
S_H = 295.08cal/mol/K
S_L = -19.59cal/mol/K
ψ_EA =12.254cal/K/mol

ΔH_A = Param(11822.57; units=cal/mol, bounds=(1e5, 2e5))
ΔH_H = Param(90678.50; units=cal/mol, bounds=(0.5e5, 1e5))
ΔH_L = Param(-6603.414; units=cal/mol, bounds=(0.5e5, 1e5))
# R = 1.987cal/K/mol
# p = K(25.0°C)*exp((ψ_EA - ΔH_A/K(25.0°C))/R) # issues
# T_ref = K(25.0°C)
p27 = Param(0.75; bounds=(0.0, 1.0)) # this is higher than 0.17 as inactivation in more in Barfield pars
T_ref = K(27.0°C)
# Have to use `withunits` here to extract the value, as Unitful and AbstractNumbers 
# don't play well together. Probably needs a fix in AbstractNumbers
Thalf_L = Param(ustrip(K, withunits(ΔH_L)/S_L); units=K, bounds=(290, 350))
Thalf_H = Param(ustrip(K, withunits(ΔH_H)/S_H); units=K, bounds=(290, 350))

growthmodel = SchoolfieldIntrinsicGrowth(p27, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)
growthresponse = Layer(:surface_temp, K, growthmodel)

GrowthMaps.rate(stripparams(growthresponse), T_ref)

#' If these are only estimated parameters, We can fit the model to a a dataset of growth
#' rate and temperature. First extract our independent and dependent variables from the example CSV:

# using CSV, DataFrames, Tables
# mkpath(joinpath(basedir, "build"))
# obsdata = CSV.File(joinpath(basedir, "swd_ecophys_data.csv"), select=[:x_value, :y_value, :y_key]) |> DataFrame
# obsdata = filter(d -> d.y_key == "r_m", obsdata) |> dropmissing

#' Then extract the required data colummns, and convert temperature values from
#' unitless Celcius to explit Kelvins, using Unitful.jl:

# obsrate = obsdata.y_value
# obstemp = obsdata.x_value * °C .|> K
# obs = collect(zip(obstemp, obsrate))

#' Now plot the growth rate curve:

temprangeC = collect(1.0:1.0:40.0)°C
temprange = K.(temprangeC)
# GrowthMaps.rate.(Ref(growthresponse), temprange)
dev = (x -> GrowthMaps.rate(stripparams(growthresponse), x)).(temprange)
p = plot(temprangeC, dev; label=false, ylabel="growth rate (1/d)", xlab = "temperature")
# scatter!(p, obs; label="observed ")
savefig(joinpath(basedir, "output/development_response.png"))

#' We can try tweaking the fitting the model manually in a user interface.
#' Model components are immutable (for performance reasons), so we wrap the model
#' in a mutable wraper so we can use the results.
#' We parametrise the model over the same temperature range that we are plotting,
#' using the :surface_temp key that the model requires:

model = Model(growthresponse)
tempdata=(surface_temp=temprange,)
interface = manualfit!(model, tempdata) #; obs=obs)

# Make an electron window:
using Blink
w = Blink.Window()
body!(w, interface)

#' If you are happy with the result, you we can update extract the manual fit
#' to use to generate our growth rate maps:

# fittedgrowth = wrapper.model

#' Note that `manualfit!` will also work for a tuple of model components
#' that use the same source data, like `(growth, heatstress, coldstress)`.
#' 
#' 
#' ## Load spatial data
#' 
#' Using the monthly aggregated SMAP dataset

# smapfolder = "/home/raf/Data/SMAP/SMAP_L4_SM_gph_v4"
smapfolder = "/home/raf/Work/Data/SMAP"
# smapfolder = "F:/SMAP/SMAP_monthly_midpoint_subset"
# smapfolder = "E:/SMAP/SMAP_L4_SM_gph_v4"

# This does the same as the above
days = (1, 8, 16, 24)
# days = 1:31
seriesall = SMAPseries(smapfolder)
series = seriesall#[Where(t -> (t >= DateTime(year) && t < DateTime(year)) && dayofmonth(t) in days)]
 
#' We can plot a layer from a file at some date in the series:
series[2][:surface_temp] |> plot

#' # Parametrising models using interactive maps
#' 
#' To find parameter values that satisfactorily explain the distribution we must
#' first aggregate the spatial data to a more manageable size so that spatial
#' predictions are more responsive to parameter adjustments.

agg = 10
aggseries = GeoData.aggregate(Center(), series, (Lon(agg), Lat(agg)); keys=(:land_fraction_wilting, :surface_temp))
p = first(aggseries[1]) |> plot

#' You can experiment with the `agg` size to compromise between quality and render time.
#' Larger values will look pixelated but will run faster.
#' 
#' First we'll run this basic growth model, calculating the daily growth rate across the year:

growthrates = mapgrowth(growthresponse;
    series=aggseries,
    tspan=DateTime(year, 1):Year(1):DateTime(year, 12)
)
plot(growthrates[Ti(1)], clim=(0, 0.15))

#' As the population growth response is only based on temperature, the annual
#' growth rate is estimated to be positive even in stressful dry environements t
#' hat have low moisture availability. These prediction can be impoved by
#' incorporating population mortality under acute stress. Now define some stressors:

coldthresh = Param(ustrip(K, 13.0°C); units=K, bounds=(-10, 20))
coldmort = Param(-0.2; units=K^-1, bounds=(-1, 0))
coldstress = Layer(:surface_temp, K, LowerStress(coldthresh, coldmort))

heatthresh = 40.0°C |> K
heatmort = Param(-0.02; units=K^-1, bounds=(-1, 0))
heatstress = Layer(:surface_temp, K, UpperStress(heatthresh, heatmort))

wiltthresh = Param(.5; bounds=(0, 1))
wiltmort = Param(-0.2; bounds=(-1, 0))
wiltstress = Layer(:land_fraction_wilting, UpperStress(wiltthresh, wiltmort));

#' check that yemen can be fit

yemen = Lon(Between(-30.0, 55.0)), Lat(Between(-30, 20))
stress = mapgrowth(growthresponse, coldstress, wiltstress, heatstress;
    series=aggseries,
    tspan=DateTime(year, 1):Year(1):DateTime(year, 12)
)
plot(mean(stress[yemen...], dims=Ti), clim=(0, 0.15))

#' To build a more complex model, we can chain components together in a tuple,
#' and run them:

tspan = DateTime(year, 1):Year(1):DateTime(year, 12)
length(tspan)
annualgrowthrates = mapgrowth((growthresponse, wiltstress, heatstress, coldstress);
    series=aggseries,
    tspan=tspan,
)
plot(annualgrowthrates[Ti(1)], clim=(0, 0.15))

#' We can also get monthly growth rates for a spatial subsection defined by longitude and latitudes:

tspan = DateTime(year):Month(1):DateTime(year, 12)
length(tspan)
monthlygrowthrates = mapgrowth((growthresponse, wiltstress, heatstress, coldstress);
    series=aggseries, tspan=tspan,
)
aust = Lon(Between(113, 153)), Lat(Between(-43, -11))
plot(monthlygrowthrates[Ti(1:3:12), aust...], clim=(0, 0.15))

#' We can also take the mean of these monthly layers and plot them:

plot(mean(monthlygrowthrates, dims=Ti), clim=(0, 0.15))

#' This looks similar to the previous yearly simulation and, indeed,
#' is equivalent due to the additive nature of exponents (e.g. exp(a) * exp(b) = exp(a + b)):

growthdiff = annualgrowthrates[Ti(1)] .- mean(monthlygrowthrates, dims=Ti)[1]
sum(map(x -> isnan(x) ? zero(x) : x, growthdiff))
tspan = DateTime(year):Month(1):DateTime(year, 12)
length(tspan)
growthrates = mapgrowth((growthresponse, heatstress, coldstress, wiltstress);
    series=aggseries, tspan=tspan,
)
plot(growthrates[Ti(1)], clims=(0.0, 0.15), axis=false)


#' ## Compare with observation data
#' 
#' To compare out simulation with observations data, we'll load them
#' from a CSV file:

csvfilename = joinpath(basedir, "data/wang2020.csv")
obs = CSV.File(csvfilename)
occurrence = collect(zip(obs.longitude, obs.latitude))

#' And scatter them on the growthrates map:

p = plot(mean(growthrates, dims=Ti), clims=(0, 0.15))
scatter!(p, occurrence;markersize=2.0, markercolor=:cyan, markershape=:cross, label=false)
display(p)
savefig(joinpath(basedir, "output/population_occurrence.png"))

#' 
#' Load monthly distibution data and plot.

# northam = Lon(Between(-125.0, -30)), Lat(Between(-50, 50))
# monthlyocc = CSV.File("data/americas_monthly.csv")
# xs = Vector(monthlyocc.longitude)
# ys = Vector(monthlyocc.latitude)
# months = Vector(monthlyocc.month)

# p = plot(monthlygrowthrates[Ti(1:12), northam...]; size=(1500, 1000), clims=(0, 0.15));
# for m in unique(months)
#     mx = xs[months .== m]
#     my = ys[months .== m]
#     scatter!(p[m], mx, my, markersize=2.0, markercolor=:cyan, markershape=:cross, label=false)
# end
# display(p)
# savefig("output/population_occurrence_monthly.png")

#' 
#' 
#' Fit a model across the whole year to see where population growth is on average positive across the year.

# growthrates = mapgrowth((growthresponse, coldstress, heatstress, wiltstress);
#     series=aggseries,
#     tspan = DateTime(2019):Year(1):DateTime(year, 12)
# )
# p = plot(growthrates[Ti(1)], clims=(0.0, 0.15))
# scatter!(p, occurrence;
#     markersize=2.0,
#     markercolor=:cyan,
#     markershape=:cross,
#     markeropacity=0.5,
#     label="obs",
# )

#' 
#' We can also look at seasons:

growthrates = mapgrowth((growthresponse, coldstress, wiltstress);
    series=aggseries,
    tspan=DateTime(year):Month(1):DateTime(year, 12),
)
p = plot(growthrates[Ti([1,4,7,10])], clims=(0.0, 0.15), axis=false)

#' 
#' As the combination of model is additive, we can prebuild parts of the model
#' we don't want to fit manually, which simplifies the interfaces and helps performance.
#' Seeing we allready fit the growth response to empiracle data, lets just fit
#' the stress responses to the map:

modelkwargs = (
    series=aggseries,
    tspan=DateTime(year):Year(1):DateTime(year, 12),
)
precomputed = mapgrowth(growthresponse; modelkwargs...)

#' Then fit the other components. `throttle` will dictate how fast the interface updates. 
#' Make it larger on a slow machine, smaller on a faster one.

# usa = Lon(Between(-125.0, -66.96)), Lat(Between(20.0, 50))

# model = Model((wiltstress, coldstress, heatstress))
# throttle = 0.2
# interface = mapfit!(model, modelkwargs;
#     window=(Ti(1), aust...),
#     # occurrence=occurrence,
#     precomputed=precomputed,
#     throttle=throttle,
#     markersize=1.0,
#     markershape=:cross,
#     markercolor=:cyan,
#     markeropacity=0.4
# )
# # Make an electron window:
# using Blink
# w = Blink.Window()
# body!(w, interface)

# And get the updated model components from the wrapper:
# wiltstress, coldstress, heatstress = parent(model)

#' Now we will put together decent population growth maps for use in other simulations using the higher resolution data, with a monthly timestep:
#' 
#' Fit a model across the whole year to see where population growth is on average positive across the year. Let's save it as a NetCDF file:

series = seriesall[Where(t -> (t >= DateTime(year) && t < DateTime(year + 1)) && dayofmonth(t) in days)]
growthrates = mapgrowth((growthresponse, wiltstress, coldstress, heatstress),
    series=series,
    tspan=DateTime(year, 1):Month(1):DateTime(year, 12)
);
for t in eachindex(dims(growthrates, Ti))
    path = joinpath(basedir, "output/growthrates")
    mkpath(path)
    tpad = lpad(t, 2, '0')
    write(joinpath(path, "growthrates_$tpad.tif"), GDALarray, growthrates[Ti(t)])
end

#' This is ready to use for estimating establishment potential, or in a dispersal
#' simulation using Dispersal.jl.
