#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using PhoenixVoltages.Fitting
using MAT
using LinearAlgebra

const centers = Solutions.CenterTracker()
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)

const wide_size = 37
const medium_size = 15
const narrow_size = 5

const fitter_wide = Fitting.PolyFitter(2, 2, 8, sizes=(5, 5, wide_size * 2 + 1))
const fitter_medium = Fitting.PolyFitter(2, 2, 6, sizes=(5, 5, medium_size * 2 + 1))
const fitter_narrow = Fitting.PolyFitter(2, 2, 4, sizes=(5, 5, narrow_size * 2 + 1))
const fitter_local = Fitting.PolyFitter(2, 2, 2)

const fit_cache_wide = Potentials.FitCache(fitter_wide, solution)
const fit_cache_medium = Potentials.FitCache(fitter_medium, solution)
const fit_cache_narrow = Potentials.FitCache(fitter_narrow, solution)
const fit_cache_local = Potentials.FitCache(fitter_local, solution)

# Plan
# [-3045, -250]: moving fitter_medium @ 3x+0 offset, center fitter_wide @ 1x+0 offset
# [-249, -120]: moving fitter_narrow @ 3x+0 offset, center fitter_medium @ 1x+0 offset
# [-119, -80]: moving fitter_narrow @ 3x+offset, center fitter_medium @ 1x+offset
#              offset linearly increasing until fitting range at the ion position
#              by the end of the range
# [-79, -35]: moving fitter_narrow @ edge, center fitter_medium @ 1x+edge
#             Trap potential for moving one linearly relaxing til 1x
# [-34, -5]: moving fitter_narrow @ 1x+edge, center fitter_medium @ 1x+edge
# [-4, -1]: moving fitter_medium @ 1x+edge, center fitter_medium @ 1x+edge
# [0, 0]: single center fitter_wide @ 1x+0offset
#
# From the second group onward, we will also start to care about the center voltage.
# In the second group, we'll gradually limit the amount of voltage difference we allow.
# From the third group on, we'll require the voltages to be identical.

const loading_um = -3045
const simple_end_um = -250
const simple2_end_um = -120
const center_shift_end_um = -80
const relax_end_um = -35
const relaxed_move_end_um = -5
const relaxed_move2_end_um = -1
const center_um = 0

const move_x2_init = 3
const center_x2 = 1

const stride_um = solution.stride .* 1000

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

const center_pos_idx = get_rf_center(center_um)
const center_pos_idx_r = (center_pos_idx[3], center_pos_idx[2], center_pos_idx[1])
const center_electrodes = Mappings.find_electrodes(solution.electrode_index,
                                                   center_um, min_num=20,
                                                   min_dist=350)

const prefix = joinpath(@__DIR__, "../data/merge_coeff_20220929")

function get_electrodes(xpos_um)
    eles = Mappings.find_electrodes(solution.electrode_index,
                                    xpos_um, min_num=20, min_dist=350)
    union!(eles, center_electrodes)
    return sort!(collect(eles))
end

function get_merge_coeff_data_simple(xpos_um)
    eles = get_electrodes(xpos_um)
    pos = get_rf_center(xpos_um)
    pos_r = (pos[3], pos[2], pos[1])
    neles = length(eles)
    coeff = Matrix{Float64}(undef, 19, neles)
    for i in 1:neles
        ele = eles[i]
        fit_move = get(fit_cache_medium, ele, pos_r)
        fit_center = get(fit_cache_wide, ele, center_pos_idx_r)
        coeff[1:9, i] .= Tuple(Solutions.get_compensate_terms1_nozx(fit_move, stride_um))
        coeff[10:19, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
    end

    target = zeros(size(coeff, 1))
    target[7] = move_x2_init # move-x2
    target[9 + 8] = center_x2 # center-x2
    x0 = coeff \ target
    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data_simple2(xpos_um)
    eles = get_electrodes(xpos_um)
    pos = get_rf_center(xpos_um)
    pos_r = (pos[3], pos[2], pos[1])
    neles = length(eles)
    coeff = Matrix{Float64}(undef, 20, neles)
    for i in 1:neles
        ele = eles[i]
        fit_move = get(fit_cache_narrow, ele, pos_r)
        fit_center = get(fit_cache_medium, ele, center_pos_idx_r)

        move_0 = get(fit_cache_local, ele, pos_r)[0, 0, 0]
        center_0 = get(fit_cache_local, ele, center_pos_idx_r)[0, 0, 0]

        coeff[1:9, i] .= Tuple(Solutions.get_compensate_terms1_nozx(fit_move, stride_um))
        coeff[10:19, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
        coeff[20, i] = move_0 - center_0
    end

    target = zeros(size(coeff, 1))
    target[7] = move_x2_init # move-x2
    target[9 + 8] = center_x2 # center-x2
    x0 = coeff \ target

    fill!(target, 0)
    target[20] = (simple2_end_um - xpos_um + 1) / (simple2_end_um - simple_end_um)
    x_offset = coeff \ target

    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "limited_solution"=>x_offset,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data_center_shift(xpos_um)
    eles = get_electrodes(xpos_um)
    pos = get_rf_center(xpos_um)
    pos_r = (pos[3], pos[2], pos[1])

    fit_shift_scale = (xpos_um - simple2_end_um) / (center_shift_end_um - simple2_end_um)
    fit_pos_r = pos_r .- (0, 0, narrow_size * fit_shift_scale)
    fit_center_pos_r = center_pos_idx_r .+ (0, 0, medium_size * fit_shift_scale)

    neles = length(eles)
    coeff = Matrix{Float64}(undef, 20, neles)
    for i in 1:neles
        ele = eles[i]
        fit_move = get(fit_cache_narrow, ele, pos_r, fit_center=fit_pos_r)
        fit_center = get(fit_cache_medium, ele, center_pos_idx_r,
                         fit_center=fit_center_pos_r)

        move_0 = get(fit_cache_local, ele, pos_r)[0, 0, 0]
        center_0 = get(fit_cache_local, ele, center_pos_idx_r)[0, 0, 0]

        coeff[1:9, i] .= Tuple(Solutions.get_compensate_terms1_nozx(fit_move, stride_um))
        coeff[10:19, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
        coeff[20, i] = move_0 - center_0
    end

    target = zeros(size(coeff, 1))
    target[7] = move_x2_init # move-x2
    target[9 + 8] = center_x2 # center-x2
    x0 = coeff \ target

    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data_relax(xpos_um)
    eles = get_electrodes(xpos_um)
    pos = get_rf_center(xpos_um)
    pos_r = (pos[3], pos[2], pos[1])

    fit_pos_r = pos_r .- (0, 0, narrow_size)
    fit_center_pos_r = center_pos_idx_r .+ (0, 0, medium_size)

    neles = length(eles)
    coeff = Matrix{Float64}(undef, 20, neles)
    for i in 1:neles
        ele = eles[i]
        fit_move = get(fit_cache_narrow, ele, pos_r, fit_center=fit_pos_r)
        fit_center = get(fit_cache_medium, ele, center_pos_idx_r,
                         fit_center=fit_center_pos_r)

        move_0 = get(fit_cache_local, ele, pos_r)[0, 0, 0]
        center_0 = get(fit_cache_local, ele, center_pos_idx_r)[0, 0, 0]

        coeff[1:9, i] .= Tuple(Solutions.get_compensate_terms1_nozx(fit_move, stride_um))
        coeff[10:19, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
        coeff[20, i] = move_0 - center_0
    end

    relax_scale = (xpos_um - center_shift_end_um) / (relax_end_um - center_shift_end_um)
    target = zeros(size(coeff, 1))
    target[7] = move_x2_init + (center_x2 - move_x2_init) * relax_scale # move-x2
    target[9 + 8] = center_x2 # center-x2
    x0 = coeff \ target

    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data_relaxed_move(xpos_um)
    eles = get_electrodes(xpos_um)
    pos = get_rf_center(xpos_um)
    pos_r = (pos[3], pos[2], pos[1])

    fit_pos_r = pos_r .- (0, 0, narrow_size)
    fit_center_pos_r = center_pos_idx_r .+ (0, 0, medium_size)

    neles = length(eles)
    coeff = Matrix{Float64}(undef, 20, neles)
    for i in 1:neles
        ele = eles[i]
        fit_move = get(fit_cache_narrow, ele, pos_r, fit_center=fit_pos_r)
        fit_center = get(fit_cache_medium, ele, center_pos_idx_r,
                         fit_center=fit_center_pos_r)

        move_0 = get(fit_cache_local, ele, pos_r)[0, 0, 0]
        center_0 = get(fit_cache_local, ele, center_pos_idx_r)[0, 0, 0]

        coeff[1:9, i] .= Tuple(Solutions.get_compensate_terms1_nozx(fit_move, stride_um))
        coeff[10:19, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
        coeff[20, i] = move_0 - center_0
    end

    target = zeros(size(coeff, 1))
    target[7] = center_x2 # move-x2
    target[9 + 8] = center_x2 # center-x2
    x0 = coeff \ target

    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data_relaxed_move2(xpos_um)
    eles = get_electrodes(xpos_um)
    pos = get_rf_center(xpos_um)
    pos_r = (pos[3], pos[2], pos[1])

    fit_pos_r = pos_r .- (0, 0, medium_size)
    fit_center_pos_r = center_pos_idx_r .+ (0, 0, medium_size)

    neles = length(eles)
    coeff = Matrix{Float64}(undef, 20, neles)
    for i in 1:neles
        ele = eles[i]
        fit_move = get(fit_cache_medium, ele, pos_r, fit_center=fit_pos_r)
        fit_center = get(fit_cache_medium, ele, center_pos_idx_r,
                         fit_center=fit_center_pos_r)

        move_0 = get(fit_cache_local, ele, pos_r)[0, 0, 0]
        center_0 = get(fit_cache_local, ele, center_pos_idx_r)[0, 0, 0]

        coeff[1:9, i] .= Tuple(Solutions.get_compensate_terms1_nozx(fit_move, stride_um))
        coeff[10:19, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
        coeff[20, i] = move_0 - center_0
    end

    target = zeros(size(coeff, 1))
    target[7] = center_x2 # move-x2
    target[9 + 8] = center_x2 # center-x2
    x0 = coeff \ target

    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data_center(xpos_um)
    @assert xpos_um == 0
    eles = get_electrodes(xpos_um)
    neles = length(eles)
    coeff = Matrix{Float64}(undef, 10, neles)
    for i in 1:neles
        ele = eles[i]
        fit_center = get(fit_cache_wide, ele, center_pos_idx_r)
        coeff[1:10, i] .= Tuple(Solutions.get_compensate_terms1(fit_center, stride_um))
    end

    target = zeros(size(coeff, 1))
    target[8] = center_x2 # center-x2
    x0 = coeff \ target
    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]

    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

function get_merge_coeff_data(xpos_um)
    @show xpos_um
    if xpos_um <= simple_end_um
        return get_merge_coeff_data_simple(xpos_um)
    elseif xpos_um <= simple2_end_um
        return get_merge_coeff_data_simple2(xpos_um)
    elseif xpos_um <= center_shift_end_um
        return get_merge_coeff_data_center_shift(xpos_um)
    elseif xpos_um <= relax_end_um
        return get_merge_coeff_data_relax(xpos_um)
    elseif xpos_um <= relaxed_move_end_um
        return get_merge_coeff_data_relaxed_move(xpos_um)
    elseif xpos_um <= relaxed_move2_end_um
        return get_merge_coeff_data_relaxed_move2(xpos_um)
    else
        return get_merge_coeff_data_center(xpos_um)
    end
end

const xpos_ums = loading_um:center_um
const coeff_data = [@time(get_merge_coeff_data(xpos_um)) for xpos_um in xpos_ums]

matopen("$(prefix).mat", "w") do mat
    write(mat, "data", coeff_data)
    write(mat, "electrode_names", solution.electrode_names)
end
