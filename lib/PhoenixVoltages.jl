#!/usr/bin/julia

module PhoenixVoltages

function gradient end

include("mappings.jl")
include("outputs.jl")
include("fitting.jl")
include("potentials.jl")
include("optimizers.jl")
include("solutions.jl")

end
