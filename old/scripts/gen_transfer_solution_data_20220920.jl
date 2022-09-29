#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT

const coeff_data_nozx, electrode_names_nozx = matopen(joinpath(@__DIR__, "../data/compensate_nozx_20220920.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const coeff_data_zx, electrode_names_zx = matopen(joinpath(@__DIR__, "../data/compensate_20220920.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const xpos_ums = [data["xpos_um"] for data in coeff_data_zx]

@assert size(coeff_data_zx) == size(coeff_data_nozx)
@assert electrode_names_zx == electrode_names_nozx
@assert [data["xpos_um"] for data in coeff_data_nozx] == xpos_ums

const electrode_names = electrode_names_nozx

const prefix = joinpath(@__DIR__, "../data/transfer_20220920")

function solve_transfer(i)
    xpos_um = xpos_ums[i]
    @show xpos_um
    # Stitch the zx vs no-zx solutions together at around -950 um,
    # where the effect of relaxing the zx term should be minimal
    # at least on the X2 anx YZ terms.
    if xpos_um < -950
        data = coeff_data_nozx[i]
        solution = data["solution"]
        x0 = solution[:, 7] .+ solution[:, 5] .* 0.75
    else
        data = coeff_data_zx[i]
        solution = data["solution"]
        x0 = solution[:, 8] .+ solution[:, 5] .* 0.75
    end
    @assert data["xpos_um"] == xpos_um
    vals = Optimizers.optimize_minmax_span(data["free_solution"], x0)
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>xpos_um)
end
const transfer_solutions = [@time(solve_transfer(i)) for i in 1:length(xpos_ums)]

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
