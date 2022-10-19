#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
using PhoenixVoltages.Solutions
using PhoenixVoltages.Solutions: l_unit, V_unit, l_unit_um, V_unit_uV
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using PhoenixVoltages.Fitting
using LinearAlgebra
using PyPlot
using NaCsPlot

const centers = Solutions.CenterTracker()
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)

const loading_um = -3045
const center_um = 0
const stride_um = solution.stride .* 1000

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end
const center_pos_idx = get_rf_center(center_um)

const center_electrodes = sort!(collect(
    Mappings.find_electrodes(solution.electrode_index,
                             center_um, min_num=4, min_dist=0)))

function get_fits(order, xsize)
    fitter = Fitting.PolyFitter(2, 2, order, sizes=(5, 5, xsize))
    fit_cache = Potentials.FitCache(fitter, solution)

    res = Matrix{Float64}(undef, length(center_electrodes), 4)
    center_r = (center_pos_idx[3], center_pos_idx[2], center_pos_idx[1])
    for (i, ele) in enumerate(center_electrodes)
        fit = get(fit_cache, ele, center_r)
        # terms = Solutions.get_compensate_terms1(fit, stride_um)
        res[i, 1] = order >= 4 ? fit[0, 0, 4] : 0
        res[i, 2] = order >= 5 ? fit[0, 0, 5] : 0
        res[i, 3] = order >= 6 ? fit[0, 0, 6] : 0
        res[i, 4] = order >= 7 ? fit[0, 0, 7] : 0
    end
    return vec(res)
end

function get_all_fits(order, xsizes)
    res = Matrix{Float64}(undef, length(xsizes), 4 * length(center_electrodes))
    for (i, xsize) in enumerate(xsizes)
        res[i, :] .= get_fits(order, xsize)
    end
    return res
end

const xsizes = 15:400

figure(figsize=[6.4 * 3, 4.8 * 3])
for (i, order) in enumerate(4:10)
    subplot(3, 3, i)
    plot(xsizes, get_all_fits(order, xsizes))
    grid()
    title("Order: $order")
end
tight_layout()

NaCsPlot.maybe_show()
