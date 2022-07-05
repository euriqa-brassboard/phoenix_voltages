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
using LsqFit
using LinearAlgebra

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Dict{String,String}()

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = Solutions.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_compensate_terms1(xpos_um)
    @show xpos_um
    return Solutions.get_compensate_terms1(fits_cache, get_rf_center(xpos_um))
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

const xpos_ums = -1715:-1690

const electrodes = ["S" .* string.(0:11); "O0"; "O1";
                    "Q44"; "Q45"; "Q0"; "Q1"; "Q2"; "Q3"]
const electrode_ids = sort!(unique(getindex.(Ref(solution.electrode_index),
                                             electrodes)))
const comp_fits2 = [[get_compensate_fits1(xpos_um, ele) for xpos_um in xpos_ums]
                    for ele in electrode_ids]

const comp_terms2 = [[Solutions.get_compensate_terms1(
    fit, solution.stride .* 1000) for fit in fits] for fits in comp_fits2]

function normalize_1(ary)
    return ary ./ maximum(abs.(ary))
end

# figure(figsize=[6.4 * 3, 4.8 * 3])
# for (i, fld) in enumerate([:dx, :dy, :dz, :xy, :yz, :zx, :x3, :x4])
#     subplot(3, 3, i)
#     for (ele, terms) in zip(electrodes, comp_terms2)
#         plot(xpos_ums, normalize_1([getfield(term, fld) / term.x2 for term in terms]),
#              label=ele)
#     end
#     grid()
#     axvline(-1702)
#     legend(ncol=2, fontsize=8)
#     title("$fld")
# end

# figure(figsize=[6.4 * 3, 4.8 * 3])
# for (i, fld) in enumerate([:x2, :dx, :dy, :dz, :xy, :yz, :zx, :x3, :x4])
#     subplot(3, 3, i)
#     for (ele, terms) in zip(electrodes, comp_terms2)
#         plot(xpos_ums, normalize_1([getfield(term, fld) for term in terms]),
#              label=ele)
#     end
#     grid()
#     axvline(-1702)
#     legend(ncol=2, fontsize=8)
#     title("$fld")
# end

figure(figsize=[6.4 * 2, 4.8 * 1])
for (i, fld) in enumerate([:dx, :zx])
    subplot(1, 2, i)
    for (ele, terms) in zip(electrodes, comp_terms2)
        plot(xpos_ums, normalize_1([getfield(term, fld) / term.x2 for term in terms]),
             label=ele)
    end
    grid()
    axvline(-1702)
    legend(ncol=2, fontsize=8)
    title("$fld")
end

figure(figsize=[6.4 * 2, 4.8 * 1])
for (i, fld) in enumerate([:dx, :zx])
    subplot(1, 2, i)
    for (ele, terms) in zip(electrodes, comp_terms2)
        plot(xpos_ums, normalize_1([getfield(term, fld) for term in terms]),
             label=ele)
    end
    grid()
    axvline(-1702)
    legend(ncol=2, fontsize=8)
    title("$fld")
end

NaCsPlot.maybe_show()
