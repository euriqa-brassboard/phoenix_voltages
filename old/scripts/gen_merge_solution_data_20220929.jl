#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/merge_coeff_20220929.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const xpos_ums = [data["xpos_um"] for data in coeff_data]

const prefix = joinpath(@__DIR__, "../data/merge_20220929")

function solve_transfer(i)
    xpos_um = xpos_ums[i]
    @show xpos_um
    data = coeff_data[i]
    x0 = data["solution"]
    @assert data["xpos_um"] == xpos_um
    if "limited_solution" in keys(data)
        limited = data["limited_solution"]
        limited = reshape(limited, size(limited, 1), :)
        vals = Optimizers.optimize_minmax_span(data["free_solution"], x0,
                                               limited=limited)
    else
        vals = Optimizers.optimize_minmax_span(data["free_solution"], x0)
    end
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>xpos_um)
end
const transfer_solutions = [@time(solve_transfer(i)) for i in 1:length(xpos_ums)]

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
