#!/usr/bin/julia

module PhoenixVoltages

function gradient end

include("mappings.jl")
include("output_files.jl")
include("voltage_solution.jl")
include("fit_multipole.jl")
include("optimizers.jl")
include("solution_process.jl")

end
