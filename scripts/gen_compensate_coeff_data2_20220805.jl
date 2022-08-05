#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using MAT
using LinearAlgebra

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = Solutions.compensate_fitter1_2(solution)

const prefix = joinpath(@__DIR__, "../data/compensate2_20220805")

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_compensate_coeff_data(xpos_um)
    @show xpos_um
    eles, coeff = Solutions.get_compensate_coeff2(fits_cache, get_rf_center(xpos_um))
    x0 = coeff \ Matrix(I, 11, 11)
    B = qr(coeff').Q[:, (11 + 1):end]
    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

const xpos_ums = -3220:1330
const coeff_data = [@time(get_compensate_coeff_data(xpos_um)) for xpos_um in xpos_ums]

matopen("$(prefix).mat", "w") do mat
    write(mat, "data", coeff_data)
    write(mat, "electrode_names", solution.electrode_names)
end
