#!/usr/bin/julia

module PhoenixVoltages

function gradient end

include("mappings.jl")
include("outputs.jl")
include("potentials.jl")
include("fit_multipole.jl")
include("optimizers.jl")
include("solution_process.jl")

end
