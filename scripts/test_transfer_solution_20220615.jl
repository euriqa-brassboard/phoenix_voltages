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
using LsqFit
using LinearAlgebra

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

function get_compensate_terms1(xpos_um)
    @show xpos_um
    return ProcessSolution.get_compensate_terms1(fits_cache, get_rf_center(xpos_um))
end

function get_compensate_terms1(xpos_um, electrode)
    @show xpos_um
    return ProcessSolution.get_compensate_terms1(get(fits_cache, electrode,
                                                     get_rf_center(xpos_um)),
                                                 solution.stride)
end

function get_compensate_fits1(xpos_um, electrode)
    @show xpos_um
    pos = get_rf_center(xpos_um)
    return get(fits_cache, electrode, (pos[3], pos[2], pos[1]))
end

function get_all_names(ids)
    names = String[]
    for id in ids
        append!(names, solution.electrode_names[id])
    end
    sort!(names)
end

# # const xpos_ums = -3220:1330
# # const xpos_ums = -1800:-1600
const xpos_um1 = -1712
const xpos_um2 = -1702

const electrodes = ["S" .* string.(0:11); "O0"; "O1";
                    "Q44"; "Q45"; "Q0"; "Q1"; "Q2"; "Q3"]
const electrode_ids = sort!(unique(getindex.(Ref(solution.electrode_index),
                                             electrodes)))
const comp_fits1 = [get_compensate_fits1(xpos_um1, ele) for ele in electrode_ids]
const comp_fits2 = [get_compensate_fits1(xpos_um2, ele) for ele in electrode_ids]

const comp_terms_sol1 = ProcessSolution.solve_terms1(comp_fits1, solution.stride .* 1000)
const comp_terms_sol2 = ProcessSolution.solve_terms1(comp_fits2, solution.stride .* 1000)

# @show comp_terms_sol1
# @show comp_terms_sol2

@show maximum(abs.(comp_terms_sol1.x2))
@show maximum(abs.(comp_terms_sol2.x2))

# function normalize_1(ary)
#     return ary ./ maximum(abs.(ary))
# end

# for fld in [:dx, :dy, :dz, :xy, :yz, :zx, :x3, :x4]
#     figure()
#     for (ele, terms) in zip(electrodes, comp_terms2)
#         plot(xpos_ums, normalize_1([getfield(term, fld) / term.x2 for term in terms]),
#              label=ele)
#     end
#     legend()
#     title("$fld")
# end

# const comp_terms = [get_compensate_terms1(xpos_um) for xpos_um in xpos_ums]
# for term in comp_terms
#     println(join(get_all_names(term[1]), "."))
# end
# const comp_files = [ProcessSolution.compensation_to_file(solution, mapfile, term...)
#                     for term in comp_terms]

# function find_min_x2(term)
#     function shifted_vector(p)
#         return term.x2 .+ term.z2 .* p[1] .+ term.xy .* p[2] .+
#             term.yz .* p[3] .+ term.zx .* p[4]
#     end
#     function model(x, p)
#         return [maximum(abs.(shifted_vector(p)))]
#     end
#     fit = curve_fit(model, [0], [0], zeros(4))
#     return (min_x2=shifted_vector(fit.param),
#             c_z2=fit.param[1],
#             c_xy=fit.param[2],
#             c_yz=fit.param[3],
#             c_zx=fit.param[4])
# end

# const min_x2_fits = [find_min_x2(term[2]) for term in comp_terms]

# # mkpath(outputdir)

# # for i in 1:length(xpos_ums)
# #     write_file(joinpath(outputdir, "comp_$(xpos_ums[i])um.txt"), comp_files[i])
# # end

# figure()
# plot(xpos_ums, [maximum(abs.(term[2].x2)) for term in comp_terms], label="x2")
# plot(xpos_ums, [maximum(abs.(fit.min_x2)) for fit in min_x2_fits], "C0--", label="x2_min")
# plot(xpos_ums, [maximum(abs.(term[2].z2)) for term in comp_terms], label="z2")
# plot(xpos_ums, [maximum(abs.(term[2].xy)) * 2 for term in comp_terms], label="xy")
# plot(xpos_ums, [maximum(abs.(term[2].yz)) * 50 for term in comp_terms], label="yz")
# plot(xpos_ums, [maximum(abs.(term[2].zx)) / 20 for term in comp_terms], label="zx")
# plot(xpos_ums, [maximum(abs.(term[2].dx)) for term in comp_terms], label="dx")
# plot(xpos_ums, [maximum(abs.(term[2].dy)) for term in comp_terms], label="dy")
# plot(xpos_ums, [maximum(abs.(term[2].dz)) for term in comp_terms], label="dz")
# plot(xpos_ums, [maximum(abs.(term[2].x3)) / 1000 for term in comp_terms], label="x3")
# plot(xpos_ums, [maximum(abs.(term[2].x4)) / 2000 for term in comp_terms], label="x4")
# legend()

# figure()
# plot(xpos_ums, [fit.c_z2 for fit in min_x2_fits], label="z2")
# plot(xpos_ums, [fit.c_xy / 2 for fit in min_x2_fits], label="xy")
# plot(xpos_ums, [fit.c_yz / 50 for fit in min_x2_fits], label="yz")
# plot(xpos_ums, [fit.c_zx * 20 for fit in min_x2_fits], label="zx")
# legend()

NaCsPlot.maybe_show()
