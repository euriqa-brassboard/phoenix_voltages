#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.VoltageSolutions
using PhoenixVoltages.PolyFit
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
const solution = ProcessSolution.ConstraintSolution(
    VoltageSolutions.import_pillbox_64(solution_file), short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)
const outputfile = ARGS[3]

const quantum_center_xidx = ProcessSolution.x_axis_to_index(solution, 0)
const quantum_center_posidx = (quantum_center_xidx, get(centers, quantum_center_xidx)...)

# (3221.0, 4.173198055959172, 3.505595438502667)
@show quantum_center_posidx

const quantum_compensate = ProcessSolution.get_compensate_terms1(
    fits_cache, quantum_center_posidx)

# @show quantum_compensate

const compfile = ProcessSolution.compensation_to_file(solution, mapfile,
                                                      quantum_compensate...)

write_file(outputfile, compfile)
