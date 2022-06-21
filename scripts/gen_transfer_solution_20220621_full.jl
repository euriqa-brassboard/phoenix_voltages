#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.VoltageSolutions
using PhoenixVoltages.PolyFit
using PhoenixVoltages.OutputFiles
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
    VoltageSolutions.import_pillbox_64(solution_file), short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)
const outputfile = ARGS[3]

function get_rf_center(xpos_um)
    xidx = ProcessSolution.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_compensate_terms1(xpos_um)
    @show xpos_um
    return ProcessSolution.get_compensate_terms1_nozx(fits_cache, get_rf_center(xpos_um))
end

function get_transfer_line(term, s)
    data = (term[2].x2 .+ term[2].yz) .* s
    return ProcessSolution.get_data_line(solution, mapfile, term[1], data)
end

const xpos_ums = -3080:1155
const scales = fill(0.25, length(xpos_ums))
scales[1:36] = range(0.2, 0.25, 36)
const comp_terms = [get_compensate_terms1(xpos_um) for xpos_um in xpos_ums]
const lines = [get_transfer_line(term, s) for (term, s) in zip(comp_terms, scales)]
push!(lines, zeros(length(mapfile.names)))

const transfer_file = TransferFile(mapfile, lines)

write_file(outputfile, transfer_file)
