#' 
#' ## Mapping incursion-point sensitivity
#' 
#' Now we will look at plotting some likely scenarios where incursion occurs
#' at a major port, then model all possible incursion scenarios.
#' 
#' Set up the simulation:

tspan = DateTime(2020, 1):timestep:DateTime(2022, 1)
init = (population=zero(populationgrid),)
output = ArrayOutput(init;
    mask=boolmask,
    aux=(growthrates=growthrates_zeromissing,),
    tspan=tspan,
)
nreps = 10

#' 
#' Set up the output grid (which can be smaller than the sim to save time)
#' Now make some pathway plots for some key incursion points.
#' 
#' And generate the plot.

function incursionplot(output, init, missingmask, key, (lat, lon), nreps)
    println(key)
    init[:population] .= 0
    init[:population][Lat(Contains(lat)), Lon(Contains(lon))] = carrycap
    steps_established = similar(output[1][:population], Int)
    steps_established .= 0
    @set! steps_established.refdims = ()
    @set! steps_established.name = "Months established"
    for rep in 1:nreps
        sim!(output, dispersalruleset; init=init)
        println("rep: ", rep)
        map(output) do step
            steps_established .+= step[:population] .>= 1
        end
    end
    println("maximum: ", maximum(steps_established))
    plot(steps_established ./ nreps .* missingmask; color=:inferno, xlabel="", ylabel="", title=key)
end

plots = []
for (key, loc) in zip(keys(incursionpoints), incursionpoints)
    push!(plots, incursionplot(output, init, missingmask, key, loc, nreps))
end

plot(plots...; size=(1000, 600))
savefig("output/months_established.png")

#' 
#' We can also save the plots as individual figures.

for (i, key) in enumerate(keys(incursionpoints))
    plot(plots[i])
    savefig("output/months_established_$key.png")
end

#' 
#' Now loop over all locations and reps to create grid where the value of each cell
#' reflects the mean number of cells invaded for an incursion commenced at that cell.

nreps = 100
scale = 2
nreps = 1
scale = 20
cellsinvaded6 = GeoData.aggregate(Center(), populationgrid, scale)
@set! cellsinvaded6.name = "Cells invaded"
cellsinvaded12 = deepcopy(cellsinvaded6)
tspan_summer = DateTime(2020, 1), DateTime(2021, 1)
tspan_winter = DateTime(2020, 7), DateTime(2021, 7)
init = (population=zero(populationgrid),)

#' 
#' Then loop over the aggregated grid. WARNING This may take hours or a day.

tiffdir = "$(@__DIR__)/output/tiffs"
rm(tiffdir; recursive=true)
mkdir(tiffdir)
build_plots!(cellsinvaded6, cellsinvaded12, dispersalruleset, init, boolmask, missingmask, nreps, scale, timestep, tspan, season) = begin
    londim, latdim = dims(init[:population], (Lon, Lat))
    latindex = reverse(val(latdim))
    trange = tspan[1]:timestep:tspan[2]
    outputs = [ArrayOutput(init, length(trange)) for t in 1:Threads.nthreads()]
    simdata = [DynamicGrids.SimData(init, dispersalruleset, first(tspan)) for t in 1:Threads.nthreads()]
    inits = [deepcopy(init) for t in 1:Threads.nthreads()]
    for i = 1:scale:size(cellsinvaded12, 1) * scale
        acc = similar(init[:population], UInt8)
        acc .= 0
        @set! acc.missingval = 0x00
        threadinit = inits[Threads.threadid()]
        threaddata = simdata[Threads.threadid()]
        threadoutput = outputs[Threads.threadid()]
        for j = 1:scale:size(cellsinvaded12, 2) * scale
            println("i, j: ", (i, j))
            boolmask[i, j] || continue
            checkbounds(Bool, acc, i, j) || continue
            lat, lon = ArchGDAL.reproject([[londim[j], latindex[i]]], crs(init[:population]), EPSG(4326))[1]
            acc .= 0
            inits[Threads.threadid()][:population] .= 0
            inits[Threads.threadid()][:population][i, j] = carrycap
            invaded6 = 0
            invaded12 = 0
            for k in 1:nreps
                sim!(threadoutput, dispersalruleset;
                     init=threadinit, tspan=tspan, simdata=threaddata)

                invaded6 += count(x -> x > one(x), threadoutput[7][:population])
                count12 = count(x -> x > one(x), threadoutput[13][:population])
                invaded12 += count12
                acc .+= threadoutput[13][:population] .> 0
            end
            cellsinvaded6[(i - 1) ÷ scale + 1, (j - 1) ÷ scale + 1] = invaded6 / nreps
            cellsinvaded12[(i - 1) ÷ scale + 1, (j - 1) ÷ scale + 1] = invaded12 / nreps
            if invaded12 > 0
                GeoData.write("$tiffdir/cells_invaded_from_$(lon)_$(lat)_$season.tif", GDALarray, acc)
            end
        end
    end
end
bm = collect(boolmask)
@time build_plots!(cellsinvaded6, cellsinvaded12, dispersalruleset, init, bm, missingmask, nreps, scale, timestep, tspan_summer, "summer")
@time build_plots!(cellsinvaded6, cellsinvaded12, dispersalruleset, init, bm, missingmask, nreps, scale, timestep, tspan_winter, "winter")

# Replace "." in lat/lon with underscores
for filename in readdir(tiffdir)
    newfilename = replace(filename, r"(.*)\.(.*)\.(.*)\.tif" => s"\1_\2_\3.tif")
    mv(joinpath(tiffdir, filename), joinpath(tiffdir, newfilename))
end
# And save a zip for uploading to RShiny
run(`zip -r incursion_tiffs.zip $tiffdir/.`)

#' 
#' Now plot:

area12 = cellsinvaded12 .* 9*9*1e-6 .* GeoData.aggregate(missingmask, scale=scale)
plot(area12; colorbar_title="Area invaded from incursion pt. (million km²)", title="Area invaded after 12 months")
savefig("output/cellsinvaded_12months_$(nreps)reps.png")
area6 = cellsinvaded6 .* 9*9*1e-6 .* GeoData.aggregate(missingmask, scale=scale)
plot(area6; colorbar_title="Area invaded from incursion pt. (million km²)", title="Area invaded after 6 months")
savefig("output/cellsinvaded_6months_$(nreps)reps.png")

#' 
#' And save the output:

using NCDatasets
write("output/cellsinvaded_6months_$(nreps)reps.ncd", NCDarray, cellsinvaded6)
write("output/cellsinvaded_12months_$(nreps)reps.ncd", NCDarray, cellsinvaded12)

A = NCDarray("output/cellsinvaded_12months_$(nreps)reps.ncd")
A = NCDarray("output/cellsinvaded_6months_$(nreps)reps.ncd")
# plot(A)

