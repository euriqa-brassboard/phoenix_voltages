#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
# using NaCsPlot
# using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return ProcessSolution.CenterTracker(read(mat, "zy_index"))
end
const short_map = ProcessSolution.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)
const outputfile = ARGS[3]

const load_center_xidx = ProcessSolution.x_axis_to_index(solution, -3.045)
const load_center_posidx = (load_center_xidx, get(centers, load_center_xidx)...)

# (176.00000000000026, 5.535493109145054, 9.640018584312939)
# @show load_center_posidx

const load_compensate = ProcessSolution.get_compensate_terms1(
    fits_cache, load_center_posidx)

# @show load_compensate

const compfile = ProcessSolution.compensation_to_file(solution, mapfile,
                                                      load_compensate...)

write_file(outputfile, compfile)
