#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.Potentials
using PhoenixVoltages.PolyFit
# using NaCsPlot
# using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return ProcessSolution.CenterTracker(read(mat, "zy_index"))
end

const solution_file = ARGS[1]
const solution = ProcessSolution.ConstraintSolution(
    Potentials.import_pillbox_64(solution_file), Dict{String,String}())
const fits_cache = ProcessSolution.compensate_fitter1(solution)

const load_center_xidx = ProcessSolution.x_axis_to_index(solution, -3.045)
const load_center_posidx = (load_center_xidx, get(centers, load_center_xidx)...)

# (176.00000000000026, 5.535493109145054, 9.640018584312939)
@show load_center_posidx

function get_x2_fit(pos)
    pos = (pos[3], pos[2], pos[1])
    fit_l0 = get(fits_cache, "L0", pos)
    fit_l1 = get(fits_cache, "L1", pos)
    # fit_l2 = get(fits_cache, "L2", pos)
    # fit_l3 = get(fits_cache, "L3", pos)
    fit_l4 = get(fits_cache, "L4", pos)
    fit_l5 = get(fits_cache, "L5", pos)
    # fit_l6 = get(fits_cache, "L6", pos)
    # fit_l7 = get(fits_cache, "L7", pos)
    fit_l8 = get(fits_cache, "L8", pos)
    fit_l9 = get(fits_cache, "L9", pos)
    fit_o0 = get(fits_cache, "O0", pos)
    fit_o1 = get(fits_cache, "O1", pos)
    return ((fit_o0 + fit_o1) * -1.89836 + (fit_l4 + fit_l5) * -0.983939
            + (fit_l0 + fit_l1 + fit_l8 + fit_l9) * 11.1849)
end

function get_yz_fit(pos)
    pos = (pos[3], pos[2], pos[1])
    fit_l0 = get(fits_cache, "L0", pos)
    fit_l1 = get(fits_cache, "L1", pos)
    # fit_l2 = get(fits_cache, "L2", pos)
    # fit_l3 = get(fits_cache, "L3", pos)
    fit_l4 = get(fits_cache, "L4", pos)
    fit_l5 = get(fits_cache, "L5", pos)
    # fit_l6 = get(fits_cache, "L6", pos)
    # fit_l7 = get(fits_cache, "L7", pos)
    fit_l8 = get(fits_cache, "L8", pos)
    fit_l9 = get(fits_cache, "L9", pos)
    fit_o0 = get(fits_cache, "O0", pos)
    fit_o1 = get(fits_cache, "O1", pos)
    return ((fit_o0 - fit_o1) * -0.32994 + (fit_l4 - fit_l5) * 0.384836
            + (fit_l0 - fit_l1 + fit_l8 - fit_l9) * 1.07679)
end

function get_dz_fit(pos)
    pos = (pos[3], pos[2], pos[1])
    fit_l0 = get(fits_cache, "L0", pos)
    fit_l1 = get(fits_cache, "L1", pos)
    # fit_l2 = get(fits_cache, "L2", pos)
    # fit_l3 = get(fits_cache, "L3", pos)
    fit_l4 = get(fits_cache, "L4", pos)
    fit_l5 = get(fits_cache, "L5", pos)
    # fit_l6 = get(fits_cache, "L6", pos)
    # fit_l7 = get(fits_cache, "L7", pos)
    fit_l8 = get(fits_cache, "L8", pos)
    fit_l9 = get(fits_cache, "L9", pos)
    fit_o0 = get(fits_cache, "O0", pos)
    fit_o1 = get(fits_cache, "O1", pos)
    return ((fit_o0 + fit_o1) * -0.777)
end

# (dx = -3.6002460743514457, dy = -19.118230859417793, dz = -8.677483111832885,
#  xy = 0.0002822932758507978, yz = -0.0017478672925034593,
#  zx = -0.0017446592745199176, z2 = -0.1868639707506263, x2 = 0.6638776755488908,
#  x3 = 8.405660805439246e-6, x4 = 0.00010237616852424318)
const fit_x2 = get_x2_fit(load_center_posidx)
@show ProcessSolution.get_compensate_terms1(fit_x2, solution.stride .* 1000)

# (dx = -0.10318868486773158, dy = 0.799257254256326, dz = -10.660438957532131,
#  xy = -4.20022908654545e-6, yz = -0.19738174779634746,
#  zx = 0.00034145534549278494, z2 = 0.003047748416465513, x2 = 0.0004781770925519853,
#  x3 = 1.6889317230326049e-6, x4 = -1.5760966430845064e-5)
const fit_yz = get_yz_fit(load_center_posidx)
@show ProcessSolution.get_compensate_terms1(fit_yz, solution.stride .* 1000)

# (dx = -3.420749586119952, dy = -6.044868085757676, dz = -999.8853408917083,
#  xy = 3.0195747514106137e-5, yz = -0.00027013902460357273,
#  zx = -0.0005365010579121651, z2 = 0.06070376179195097, x2 = 7.924978702025439e-5,
#  x3 = 2.2192028883793418e-7, x4 = 4.140523041832714e-5)
const fit_dz = get_dz_fit(load_center_posidx)
@show ProcessSolution.get_compensate_terms1(fit_dz, solution.stride .* 1000)
