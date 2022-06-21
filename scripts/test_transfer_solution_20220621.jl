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
const short_map = Dict{String,String}()

const solution_file = ARGS[1]
const solution = ProcessSolution.ConstraintSolution(
    VoltageSolutions.import_pillbox_64(solution_file), short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)

function get_rf_center(xpos_um)
    xidx = ProcessSolution.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_compensate_terms1_nozx(xpos_um)
    @show xpos_um
    return ProcessSolution.get_compensate_terms1_nozx(fits_cache, get_rf_center(xpos_um))
end

const xpos_ums = -3220:1330
const comp_terms = [get_compensate_terms1_nozx(xpos_um) for xpos_um in xpos_ums]

figure()
plot(xpos_ums, [maximum(abs.(term[2].dx)) * 1000 for term in comp_terms], label="dx")
plot(xpos_ums, [maximum(abs.(term[2].dy)) * 1000 for term in comp_terms], label="dy")
plot(xpos_ums, [maximum(abs.(term[2].dz)) * 1000 for term in comp_terms], label="dz")
plot(xpos_ums, [maximum(abs.(term[2].x2)) for term in comp_terms], label="x2")
plot(xpos_ums, [maximum(abs.(term[2].z2)) for term in comp_terms], label="z2")
plot(xpos_ums, [maximum(abs.(term[2].xy)) for term in comp_terms], label="xy")
plot(xpos_ums, [maximum(abs.(term[2].yz)) for term in comp_terms], label="yz")
for x in -3220:70:1330
    axvline(x)
end
legend()
grid()

figure()
plot(xpos_ums, [sum(abs2.(term[2].dx)) * 1000 for term in comp_terms], label="dx")
plot(xpos_ums, [sum(abs2.(term[2].dy)) * 1000 for term in comp_terms], label="dy")
plot(xpos_ums, [sum(abs2.(term[2].dz)) * 1000 for term in comp_terms], label="dz")
plot(xpos_ums, [sum(abs2.(term[2].x2)) for term in comp_terms], label="x2")
plot(xpos_ums, [sum(abs2.(term[2].z2)) for term in comp_terms], label="z2")
plot(xpos_ums, [sum(abs2.(term[2].xy)) for term in comp_terms], label="xy")
plot(xpos_ums, [sum(abs2.(term[2].yz)) for term in comp_terms], label="yz")
for x in -3220:70:1330
    axvline(x)
end
legend()
grid()

NaCsPlot.maybe_show()
