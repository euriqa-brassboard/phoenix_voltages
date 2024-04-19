#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
import PhoenixVoltages.Fitting
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using MAT
using LinearAlgebra

using JuMP
using Ipopt

using NaCsPlot
using PyPlot

# Loading input data
const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202310.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const center_fit_cache = Solutions.compensate_fitter3(solution)

# Utility functions
function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function find_index_range(center, radius, sz)
    lb = floor(Int, center - radius)
    ub = ceil(Int, center + radius)
    lb = max(1, lb)
    ub = min(sz, ub)
    return lb:ub
end

struct ElectrodeData
    slice_data::Vector{Float64}
    center_data::Vector{Float64}
end

function ElectrodeData(fit_cache, ele, xindex_range, center_index)
    solution = fit_cache.solution
    data = solution.data
    stride_ums_zyx = solution.stride .* 1000
    data = @view(data[:, :, xindex_range, ele])
    len = length(xindex_range)
    # 1, y, z, y^2, yz, z^2
    slice_data = Matrix{Float64}(undef, 6, len)
    fitter = Fitting.PolyFitter(4, 4, sizes=(7, 7))
    scales_zyx = Solutions.l_unit_um ./ stride_ums_zyx
    for i in 1:len
        data_zy = @view(data[:, :, i])
        cache_zy = Fitting.PolyFitCache(fitter, data_zy)
        fit_zy = get(cache_zy, (center_index[3], center_index[2]))
        slice_data[1, i] = fit_zy[0, 0] / Solutions.V_unit
        slice_data[2, i] = fit_zy[0, 1] / Solutions.V_unit * scales_zyx[2]
        slice_data[3, i] = fit_zy[1, 0] / Solutions.V_unit * scales_zyx[1]
        slice_data[4, i] = fit_zy[0, 2] / Solutions.V_unit * scales_zyx[2]^2 * 2
        slice_data[5, i] = fit_zy[1, 1] / Solutions.V_unit * scales_zyx[2] * scales_zyx[1]
        slice_data[6, i] = fit_zy[2, 0] / Solutions.V_unit * scales_zyx[1]^2 * 2
    end

    fit = get(fit_cache, ele, (center_index[3], center_index[2], center_index[1]))
    terms = Solutions.get_compensate_terms2(fit, stride_ums_zyx)

    return ElectrodeData(vec(slice_data),
                         [terms.dx, terms.dy, terms.dz, terms.xy, terms.yz, terms.zx,
                          terms.z2, terms.x2, terms.x3, terms.x4,
                          fit[0, 1, 2], # x^2y
                          fit[1, 0, 2], # x^2z
                          ])
end

const slice_order_mapping = Dict((0, 0) => 1, (0, 1) => 2, (1, 0) => 3,
                                 (0, 2) => 4, (1, 1) => 5, (2, 0) => 6)
function generate_term(center_index, xindex_range, stride_x_um, order)
    @assert all(0 .<= order)
    @assert order[3] + order[2] <= 2
    param_idx = slice_order_mapping[(order[3], order[2])]
    len = length(xindex_range)
    res = zeros(6, len)
    xcenter = center_index[1]
    xscale = stride_x_um / Solutions.l_unit_um
    coeff = 1 / factorial(order[1])
    for i in 1:len
        xpos = xindex_range[i] - xcenter
        res[param_idx, i] = (xpos * xscale)^order[1] * coeff
    end
    return vec(res)
end

struct AllTermsData
    ele_data::Vector{ElectrodeData}
    center_x0::Matrix{Float64}
    center_B::Matrix{Float64}
    slice_coeff::Matrix{Float64}
    slice_terms::Matrix{Float64}
end

const xregion_radius = 65
function AllTermsData(fit_cache, center_pos_um)
    solution = fit_cache.solution
    center_index = get_rf_center(center_pos_um)
    xindex_range = find_index_range(center_index[1], xregion_radius, solution.nx)
    eles = sort!(collect(Mappings.find_electrodes(solution.electrode_index,
                                                  center_pos_um,
                                                  min_num=30, min_dist=600)))

    ele_data = [ElectrodeData(fit_cache, ele, xindex_range, center_index)
                for ele in eles]

    ncenter_coeff = @show length(ele_data[1].center_data)
    center_coeff = Matrix{Float64}(undef, ncenter_coeff, length(ele_data))
    for i in 1:length(ele_data)
        center_coeff[:, i] .= ele_data[i].center_data
    end
    center_x0 = center_coeff \ Matrix(I, ncenter_coeff, ncenter_coeff)
    center_B = qr(center_coeff').Q[:, (ncenter_coeff + 1):end]

    stride_x_um = solution.stride[3] * 1000
    term_0 = generate_term(center_index, xindex_range, stride_x_um, (0, 0, 0))
    slice_coeff = Matrix{Float64}(undef, length(term_0), length(ele_data) + 1)
    slice_coeff[:, 1] .= term_0
    for i in 1:length(ele_data)
        slice_coeff[:, 1 + i] .= ele_data[i].slice_data
    end
    term_x1_raw = generate_term(center_index, xindex_range, stride_x_um, (1, 0, 0))
    term_x2_raw = generate_term(center_index, xindex_range, stride_x_um, (2, 0, 0))
    term_x3_raw = generate_term(center_index, xindex_range, stride_x_um, (3, 0, 0))
    term_x4_raw = generate_term(center_index, xindex_range, stride_x_um, (4, 0, 0))

    term_y1_raw = generate_term(center_index, xindex_range, stride_x_um, (0, 1, 0))
    term_y2_raw = generate_term(center_index, xindex_range, stride_x_um, (0, 2, 0))

    term_z1_raw = generate_term(center_index, xindex_range, stride_x_um, (0, 0, 1))
    term_z2_raw = generate_term(center_index, xindex_range, stride_x_um, (0, 0, 2))

    term_x1y1_raw = generate_term(center_index, xindex_range, stride_x_um, (1, 1, 0))
    term_y1z1_raw = generate_term(center_index, xindex_range, stride_x_um, (0, 1, 1))
    term_x1z1_raw = generate_term(center_index, xindex_range, stride_x_um, (1, 0, 1))

    term_x1y2_raw = generate_term(center_index, xindex_range, stride_x_um, (1, 2, 0))
    term_x2y2_raw = generate_term(center_index, xindex_range, stride_x_um, (2, 2, 0))
    term_x1z2_raw = generate_term(center_index, xindex_range, stride_x_um, (1, 0, 2))
    term_x2z2_raw = generate_term(center_index, xindex_range, stride_x_um, (2, 0, 2))

    term_dx = term_x1_raw .* (Solutions.l_unit_um / Solutions.V_unit_uV)
    term_dy = term_y1_raw .* (Solutions.l_unit_um / Solutions.V_unit_uV)
    term_dz = term_z1_raw .* (Solutions.l_unit_um / Solutions.V_unit_uV)

    term_xy = term_x1y1_raw
    term_yz = term_y1z1_raw
    term_zx = term_x1z1_raw
    term_z2 = term_z2_raw - term_y2_raw
    term_x2 = term_x2_raw .- (term_y2_raw .+ term_z2_raw) ./ 2

    term_x3 = term_x3_raw .- (term_x1y2_raw .+ term_x1z2_raw) ./ 2
    term_x4 = term_x4_raw .- (term_x2y2_raw .+ term_x2z2_raw) ./ 2

    slice_terms = hcat(term_dx, term_dy, term_dz,
                       term_xy, term_yz, term_zx, term_z2, term_x2,
                       term_x3, term_x4)
    return AllTermsData(ele_data, center_x0, center_B, slice_coeff, slice_terms)
end

function fit_term(model::Model, term_data::AllTermsData, term_idx, maxv, yz_weight;
                  relax_high_order=false)
    x0 = @view(term_data.center_x0[:, term_idx])
    B = term_data.center_B
    nx, nt = size(B)
    if true
        @variable(model, t[1:nt])
        x = @expression(model, B * t .+ x0)
        if relax_high_order
            nhigh_order = size(term_data.center_x0, 2) - 10
            @variable(model, ht[1:nhigh_order])
            x = @expression(model, x .+ term_data.center_x0[:, 11:end] * ht)
        end
        @constraint(model, -maxv .<= x .<= maxv)
    else
        x = x0
    end

    nvars = nx + 1
    @variable(model, dc_offset)
    vars = [dc_offset; x]
    target = @view(term_data.slice_terms[:, term_idx])
    err = @expression(model, term_data.slice_coeff * vars .- target)
    nvals = length(err)
    @variable(model, abserr[1:nvals])
    @constraint(model, .-abserr .<= err)
    @constraint(model, abserr .>= err)
    npoints = length(target) ÷ 6
    abserr = reshape(abserr, (6, npoints))

    @variable(model, maxerr)

    point_center = (1 + npoints) / 2
    max_offset = point_center - 1
    C = log(2) / max_offset^2
    for i in 1:npoints
        x_weight = exp(-(i - point_center)^2 * C)
        for j in 1:6
            e = abserr[j, i] * x_weight * yz_weight[j]
            @constraint(model, -maxerr <= e)
            @constraint(model, maxerr >= e)
        end
    end
    @objective(model, Min, maxerr)
    optimize!(model)
    return value.(vars)
end

# Basic input parameters
const center_pos_um = 0.0

const all_term_data = AllTermsData(center_fit_cache, center_pos_um)

const voltages_dx = fit_term(Model(Ipopt.Optimizer), all_term_data, 1, 0.003,
                               (10, 3, 3, 1, 3, 1))
const voltages_dy = fit_term(Model(Ipopt.Optimizer), all_term_data, 2, 0.003,
                               (5, 6, 3, 1, 3, 1))
const voltages_dz = fit_term(Model(Ipopt.Optimizer), all_term_data, 3, 0.003,
                               (5, 5, 10, 1, 3, 1))
const voltages_xy = fit_term(Model(Ipopt.Optimizer), all_term_data, 4, 6.0,
                               (400, 300, 6, 20, 200, 20))
const voltages_yz = fit_term(Model(Ipopt.Optimizer), all_term_data, 5, 3.2,
                             (40, 30, 60, 20, 300, 20), relax_high_order=true)
const voltages_zx = fit_term(Model(Ipopt.Optimizer), all_term_data, 6, 9.0,
                             (4, 30, 60, 40, 60, 40), relax_high_order=true)
const voltages_z2 = fit_term(Model(Ipopt.Optimizer), all_term_data, 7, 10.0,
                               (40, 6, 6, 200, 300, 2000))
const voltages_x2 = fit_term(Model(Ipopt.Optimizer), all_term_data, 8, 19.0,
                               (40, 30, 3, 300, 300, 2000))
const voltages_x3 = fit_term(Model(Ipopt.Optimizer), all_term_data, 9, 200.0,
                               (6, 150, 20, 300, 700, 200), relax_high_order=true)
const voltages_x4 = fit_term(Model(Ipopt.Optimizer), all_term_data, 10, 6000.0,
                               (4, 10, 10, 300, 300, 200), relax_high_order=true)

function get_nth_part(data, idx)
    return data[idx:6:end]
end

const term_names = ["1", "y", "z", "y2", "yz", "z2"]

@show extrema(voltages_dx[2:end])
term_dx_flat = all_term_data.slice_terms[:, 1]
true_dx = all_term_data.slice_coeff * voltages_dx
err_dx = term_dx_flat .- true_dx
@show extrema(err_dx)

@show extrema(voltages_dy[2:end])
term_dy_flat = all_term_data.slice_terms[:, 2]
true_dy = all_term_data.slice_coeff * voltages_dy
err_dy = term_dy_flat .- true_dy
@show extrema(err_dy)

@show extrema(voltages_dz[2:end])
term_dz_flat = all_term_data.slice_terms[:, 3]
true_dz = all_term_data.slice_coeff * voltages_dz
err_dz = term_dz_flat .- true_dz
@show extrema(err_dz)

@show extrema(voltages_xy[2:end])
term_xy_flat = all_term_data.slice_terms[:, 4]
true_xy = all_term_data.slice_coeff * voltages_xy
err_xy = term_xy_flat .- true_xy
@show extrema(err_xy)

@show extrema(voltages_yz[2:end])
term_yz_flat = all_term_data.slice_terms[:, 5]
true_yz = all_term_data.slice_coeff * voltages_yz
err_yz = term_yz_flat .- true_yz
@show extrema(err_yz)

@show extrema(voltages_zx[2:end])
term_zx_flat = all_term_data.slice_terms[:, 6]
true_zx = all_term_data.slice_coeff * voltages_zx
err_zx = term_zx_flat .- true_zx
@show extrema(err_zx)

@show extrema(voltages_z2[2:end])
term_z2_flat = all_term_data.slice_terms[:, 7]
true_z2 = all_term_data.slice_coeff * voltages_z2
err_z2 = term_z2_flat .- true_z2
@show extrema(err_z2)

@show extrema(voltages_x2[2:end])
term_x2_flat = all_term_data.slice_terms[:, 8]
true_x2 = all_term_data.slice_coeff * voltages_x2
err_x2 = term_x2_flat .- true_x2
@show extrema(err_x2)

@show extrema(voltages_x3[2:end])
term_x3_flat = all_term_data.slice_terms[:, 9]
true_x3 = all_term_data.slice_coeff * voltages_x3
err_x3 = term_x3_flat .- true_x3
@show extrema(err_x3)

@show extrema(voltages_x4[2:end])
term_x4_flat = all_term_data.slice_terms[:, 10]
true_x4 = all_term_data.slice_coeff * voltages_x4
err_x4 = term_x4_flat .- true_x4
@show extrema(err_x4)

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_dx_flat, i), ls="--")
    plot(get_nth_part(true_dx, i))
    grid()
    title(term_names[i])
end
suptitle("dx")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_dy_flat, i), ls="--")
    plot(get_nth_part(true_dy, i))
    grid()
    title(term_names[i])
end
suptitle("dy")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_dz_flat, i), ls="--")
    plot(get_nth_part(true_dz, i))
    grid()
    title(term_names[i])
end
suptitle("dz")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_xy_flat, i), ls="--")
    plot(get_nth_part(true_xy, i))
    grid()
    title(term_names[i])
end
suptitle("xy")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_yz_flat, i), ls="--")
    plot(get_nth_part(true_yz, i))
    grid()
    title(term_names[i])
end
suptitle("yz")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_zx_flat, i), ls="--")
    plot(get_nth_part(true_zx, i))
    grid()
    title(term_names[i])
end
suptitle("zx")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_z2_flat, i), ls="--")
    plot(get_nth_part(true_z2, i))
    grid()
    title(term_names[i])
end
suptitle("z2")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_x2_flat, i), ls="--")
    plot(get_nth_part(true_x2, i))
    grid()
    title(term_names[i])
end
suptitle("x2")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_x3_flat, i), ls="--")
    plot(get_nth_part(true_x3, i))
    grid()
    title(term_names[i])
end
suptitle("x3")

figure(figsize=[6.4 * 3, 4.8 * 2])
for i in 1:6
    subplot(2, 3, i)
    plot(get_nth_part(term_x4_flat, i), ls="--")
    plot(get_nth_part(true_x4, i))
    grid()
    title(term_names[i])
end
suptitle("x4")

NaCsPlot.maybe_show()
