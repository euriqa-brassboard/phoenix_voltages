#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.Potentials
using PhoenixVoltages.PolyFit
using NaCsPlot
using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return ProcessSolution.CenterTracker(read(mat, "zy_index"))
end
const short_map = ProcessSolution.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = ProcessSolution.ConstraintSolution(
    Potentials.import_pillbox_64(solution_file), short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const prefix = joinpath(@__DIR__, "../data/transfer_20220630")

function get_rf_center(xpos_um)
    xidx = ProcessSolution.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_transfer1(xpos_um)
    @show xpos_um
    eles, vals = ProcessSolution.get_transfer1(fits_cache, get_rf_center(xpos_um))
    ele_names = [solution.electrode_names[ele] for ele in eles]
    return [ele_names, vals]
end

const xpos_ums = -3220:1330
const transfer_solutions = [get_transfer1(xpos_um) for xpos_um in xpos_ums]

matopen("$(prefix).mat", "w") do mat
    write(mat, "xpos_ums", [xpos_ums;])
    write(mat, "transfer_solutions", transfer_solutions)
end
