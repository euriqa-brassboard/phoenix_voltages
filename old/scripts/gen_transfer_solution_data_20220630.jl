#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using NaCsPlot
using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = Solutions.compensate_fitter1(solution)

const prefix = joinpath(@__DIR__, "../data/transfer_20220630")

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function solve_transfer(xpos_um)
    @show xpos_um
    eles, vals = Solutions.solve_transfer1(fits_cache, get_rf_center(xpos_um))
    ele_names = [solution.electrode_names[ele] for ele in eles]
    return [ele_names, vals]
end

const xpos_ums = -3220:1330
const transfer_solutions = [solve_transfer(xpos_um) for xpos_um in xpos_ums]

matopen("$(prefix).mat", "w") do mat
    write(mat, "xpos_ums", [xpos_ums;])
    write(mat, "transfer_solutions", transfer_solutions)
end
