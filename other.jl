
mitchell = CSV.File("data/Mitchell 1991 - Trap counts.csv")
n = Vector(mitchell.trap_count)
t = map(x -> DateTime(x, DateFormat("d/m/Y")), Vector(mitchell.date))
t = t .+ Year(36)
minimum(t)
loc = Vector(mitchell.location)
unique(loc)
us_cline = (
    St_Hyacinth = (45.627716, -72.956093),
    Tifton_GA = (31.445379, -83.513415),
    Gainsville_FL = (29.661792, -82.326272),
    Homestead_FL=(25.470084, -80.464693),
    Puerto_Rico=(18.267668, -66.445990),
    Virgin_Islands=(17.728809, -64.821270),
    Guadeloupe=(16.187087, -61.676442),
    French_Guiana=(4.905187, -52.418892),
)
popxy = zeros(length(us_tspan))
for i in 1:length(us_cline)
    println(i)
    (y, x) = us_cline[i]
    for j in 1:length(us_tspan)
        popxy[j] = output[j][:population][Lat(Contains(y)), Lon(Contains(x))]
    end
    p=plot(us_tspan, (popxy .- minimum(popxy))./(maximum(popxy) .- minimum(popxy)),
        xticks = us_tspan[1]:Month(6):us_tspan[2], title=keys(us_cline)[i])
    iloc = keys(us_cline)[i]
    it = t[loc .== String(iloc)]
    in = n[loc .== String(iloc)]
    scatter!(p, it, in/maximum(in))
    ifilename = joinpath(basedir, "output/population_growth_$(iloc).png")
    display(p)
    savefig(ifilename)
end

mysubset=Vector{Int}(undef,12)
for c = 1:12
    println(c)
    mysubset[c] =
        findall((us_tspan .>= Date(2022,1,1) + Month(c - 1)))[1]
end
us_tspan[mysubset]

for i in mysubset
    usspreaddir = joinpath(basedir, "output/usa_spread/")
    mkpath(usspreaddir)
    fn = joinpath(usspreaddir, string(Date(us_tspan[i]), ".ncd"))
    write(fn, NCDarray, output[i][:population])
end

output = ArrayOutput((population=populationgridusa,);
    mask=boolmaskusa,
    aux=(growthrates=growthrates_zeromissingusa,),
    tspan=us_tspan,
)
@time sim!(output, dispersalruleset)
(y, x) = incursionpoints[:Tifton_GA]
plot(output[2][:population] .* missingmaskusa)

popxy = [frame[:population][Lat(Contains(y)), Lon(Contains(x))] for frame in output]
output[2][:population]
maximum(output[1][:population])
typeof(output[2][:population])
plot(us_tspan, popxy)


using DynamicGridsGtk
# init_popgrid!(populationgrid, incursionpoints[:Texas], carrycap)
output = GtkOutput(
    (population=populationgrid,),
    mask=boolmask,
    aux=(growthrates=test,),
    tspan=us_tspan,
    store=true,
    processor=inferno,
    minval=zero(carrycap),
    maxval=carrycap,
)
sim!(output, dispersalruleset)


using DataFrames, CSV, NCDatasets
# set year and index dict
iyear = "2016"
GRdict = Dict([
    ("2016", growthrates_zeromissingaus2016),
    ("2017", growthrates_zeromissingaus2017),
    ("2018", growthrates_zeromissingaus2018),
    ("2019", growthrates_zeromissingaus2019)])

for iyear in ("2016","2017","2018","2019")
    output = ArrayOutput((population=populationgridaus,);
        mask=boolmaskaus,
        aux=(growthrates=GRdict[iyear],),
        tspan=aus_tspan,
    )
    @time sim!(output, dispersalruleset)
    mysubset=Vector{Int}(undef,12)
    for c = 1:12
        println(c)
        # get monthly seasonality in third year
        mysubset[c] =
            findall((aus_tspan .>= Date(2022,1,1) + Month(c - 1)))[1]
    end
    aus_tspan[mysubset]
    for i in mysubset
        # hack date and save output
        spreaddir = joinpath(basedir, "output/aus_spread/")
        mkpath(spreaddir)
        fn = joinpath(spreaddir, string(replace(string(Date(aus_tspan[i])), "2022"=>iyear), ".ncd"))
        @show fn
        write(fn, NCDarray, output[i][:population])
    end
    GRDCregion = CSV.read(joinpath(basedir, "data/growingregion_coords.csv"))
    println(GRDCregion)
    d = DataFrame(
        GRDCregion = String[],
        rep = Int[],
        date = Date[],
        lon = Float64[],
        lat = Float64[],
        popsize = Float64[],
    )
    popxy = zeros(length(aus_tspan))
    nreps = 100
    for rep in 1:nreps
        println(rep, " of ", nreps)
        @time sim!(output, dispersalruleset)
        p=plot()
        for i in 1:size(GRDCregion, 1)
            iloc = GRDCregion.Region[i]
            x = GRDCregion.Longitude[i]
            y = GRDCregion.Latitude[i]
                for j in 1:length(aus_tspan)
                    popxy[j] = output[j][:population][Lat(Contains(y)), Lon(Contains(x))]
                    push!(d, (iloc, rep, aus_tspan[j], x, y, popxy[j]))
                end
                plot!(p, aus_tspan, (popxy .- minimum(popxy))./(maximum(popxy) .- minimum(popxy)) ,
                    xticks = aus_tspan[1]:Month(6):aus_tspan[2], label=iloc)
        end
        # display(p)
    end
    ifilename = joinpath(basedir, string("output/GRDC_region_FAW_activity_",iyear,".png"))
    savefig(ifilename)
    CSV.write(replace(ifilename,  ".png" => ".csv"), d)
end
