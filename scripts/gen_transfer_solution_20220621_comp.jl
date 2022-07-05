#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
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

const mapfile = load_file(ARGS[2], MapFile)
const outputdir = ARGS[3]

function get_rf_center(xpos_um)
    xidx = ProcessSolution.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_compensate_terms1(xpos_um)
    @show xpos_um
    return ProcessSolution.get_compensate_terms1_nozx(fits_cache, get_rf_center(xpos_um))
end

const xpos_ums = -3220:1330
const comp_terms = [get_compensate_terms1(xpos_um) for xpos_um in xpos_ums]
const comp_files = [ProcessSolution.compensation_to_file(solution, mapfile, term...)
                    for term in comp_terms]

mkpath(outputdir)

for i in 1:length(xpos_ums)
    write_file(joinpath(outputdir, "comp_$(xpos_ums[i])um.txt"), comp_files[i])
end
