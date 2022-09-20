#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/compensate_nozx_20220920.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const prefix = joinpath(@__DIR__, "../data/transfer_20220920")

function solve_transfer(data)
    @show data["xpos_um"]
    solution = data["solution"]
    x0 = solution[:, 7] .+ solution[:, 5] .* 0.75
    vals = Optimizers.optimize_minmax_span(data["free_solution"], x0)
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end
const transfer_solutions = [@time(solve_transfer(data)) for data in coeff_data]

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
