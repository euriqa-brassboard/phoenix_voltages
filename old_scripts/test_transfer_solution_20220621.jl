#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using NaCsPlot
using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))
# const short_map = Dict{String,String}()

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = Solutions.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function solve_compensate(xpos_um)
    @show xpos_um
    return Solutions.solve_compensate1_nozx(fits_cache, get_rf_center(xpos_um))
end

const xpos_ums = -3080:35:1155
const comp_terms = [solve_compensate(xpos_um) for xpos_um in xpos_ums]

const transfer_terms = [(term[2].x2 .+ term[2].yz) .* 0.25 for term in comp_terms]
transfer_terms[1] = transfer_terms[1] .* 0.8

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
plot(xpos_ums, [maximum(abs.(term)) for term in transfer_terms])
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
