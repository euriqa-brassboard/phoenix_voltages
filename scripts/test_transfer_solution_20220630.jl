#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.VoltageSolutions
using PhoenixVoltages.PolyFit
using PhoenixVoltages.Outputs
using NaCsPlot
using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return ProcessSolution.CenterTracker(read(mat, "zy_index"))
end
const short_map = ProcessSolution.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))
# const short_map = Dict{String,String}()

const solution_file = ARGS[1]
const solution = ProcessSolution.ConstraintSolution(
    VoltageSolutions.import_pillbox_64(solution_file), short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)

function get_rf_center(xpos_um)
    xidx = ProcessSolution.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_transfer1(xpos_um)
    @show xpos_um
    return ProcessSolution.get_transfer1(fits_cache, get_rf_center(xpos_um))
end

const xpos_ums = -3080:1155
const transfer_terms = [get_transfer1(xpos_um)[2] for xpos_um in xpos_ums]

figure()
for x in -3220:70:1330
    axvline(x, color="C1")
end
plot(xpos_ums, [maximum(abs.(term)) for term in transfer_terms])
grid()

NaCsPlot.maybe_show()
