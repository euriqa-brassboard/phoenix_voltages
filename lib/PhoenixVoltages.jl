#!/usr/bin/julia

module PhoenixVoltages

function gradient end

include("mappings.jl")
include("outputs.jl")
include("potentials.jl")
include("fitting.jl")
include("optimizers.jl")
include("solution_process.jl")

end
