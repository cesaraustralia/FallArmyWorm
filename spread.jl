#' 
#' # Dispersal simulations
#' In this example we will run a simulation of the spread of the Fall Armyworm _Spodoptera frugiperda_ accross the continental USA.
#' 
#' ### GTK window
#' 
#' Set up a simple rule that just copies the growthrate layer to a grid so that
#' we can view it next to the simulation.

aus_tspan = DateTime(2020, 3):timestep:DateTime(2023, 3)
populationgridaus = replace_missing(zero(growthratesaus[Ti(1)]), NaN)
# init_popgrid!(populationgrid, incursionpoints[:Seisia], 1.0e6)
populationgridaus += (growthmaskaus[Ti(1)].>0) .* carrycap

#' 
#' Run the simulation in a DynamicGridsInteract window:

using DynamicGridsInteract, DynamicGridsGtk
output = GtkOutput(
#output = ElectronOutput(
    (population=populationgridaus,);
    ruleset=dispersalruleset,
    mask=boolmaskaus,
    aux=(growthrates=growthrates_zeromissingaus,),
    tspan=aus_tspan,
    store=true,
    processor=inferno,
    minval=zero(carrycap),
    maxval=carrycap,
)
display(output)
sim!(output, dispersalruleset)

savegif(joinpath(basedir, "output/dispersalaus.gif"), output, dispersalruleset; processor=inferno, fps=3)
