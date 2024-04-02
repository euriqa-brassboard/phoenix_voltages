#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
import PhoenixVoltages.Fitting
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using MAT
using LinearAlgebra

using NaCsPlot
using PyPlot

using JuMP
using Ipopt

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202310.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)

const prefix = joinpath(@__DIR__, "../data/bestfit_20240325")

function get_rf_center_zyx(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    yidx, zidx = get(centers, xidx)
    return (zidx, yidx, xidx)
end

function find_index_range(center, radius, sz)
    lb = floor(Int, center - radius)
    ub = ceil(Int, center + radius)
    lb = max(1, lb)
    ub = min(sz, ub)
    return lb:ub
end

const center_pos_um = 0.0

const ele_select = sort!(collect(Mappings.find_electrodes(solution.electrode_index,
                                                          center_pos_um,
                                                          min_num=30, min_dist=600)))
const center_index_zyx = get_rf_center_zyx(center_pos_um)
const region_radius_zyx = (3, 3, 75)

const index_range_zyx = find_index_range.(center_index_zyx, region_radius_zyx,
                                          (solution.nz, solution.ny, solution.nx))

function get_ele_data(data, index_range_zyx, ele,
                      center_index_zyx, stride_zyx)
    data = @view(data[index_range_zyx..., ele])
    len = length(index_range_zyx[3])
    # 1, y, z, y^2, yz, z^2
    res = Matrix{Float64}(undef, 6, len)
    fitter = Fitting.PolyFitter(4, 4)
    scales_zyx = Solutions.l_unit_um ./ (stride_zyx .* 1000)
    for i in 1:len
        data_zy = @view(data[:, :, i])
        cache_zy = Fitting.PolyFitCache(fitter, data_zy)
        fit_zy = get(cache_zy, (center_index_zyx[1], center_index_zyx[2]))
        res[1, i] = fit_zy[0, 0] / Solutions.V_unit
        res[2, i] = fit_zy[0, 1] / Solutions.V_unit * scales_zyx[2]
        res[3, i] = fit_zy[1, 0] / Solutions.V_unit * scales_zyx[1]
        res[4, i] = fit_zy[0, 2] / Solutions.V_unit * scales_zyx[2]^2 * 2
        res[5, i] = fit_zy[1, 1] / Solutions.V_unit * scales_zyx[1] * scales_zyx[1]
        res[6, i] = fit_zy[2, 0] / Solutions.V_unit * scales_zyx[1]^2 * 2
    end
    return res
end

const ele_data = [get_ele_data(solution.data, index_range_zyx, ele,
                               center_index_zyx, solution.stride) for ele in ele_select]

const order_mapping = Dict((0, 0) => 1, (0, 1) => 2, (1, 0) => 3,
                           (0, 2) => 4, (1, 1) => 5, (2, 0) => 6)

function generate_term(center_index_zyx, index_range_zyx, stride_zyx, order_zyx)
    @assert all(0 .<= order_zyx)
    @assert order_zyx[1] + order_zyx[2] <= 2
    param_idx = order_mapping[(order_zyx[1], order_zyx[2])]
    len = length(index_range_zyx[3])
    res = zeros(6, len)
    xcenter = center_index_zyx[3]
    xscale = (stride_zyx[3] * 1000) / Solutions.l_unit_um
    coeff = 1 / factorial(order_zyx[3])
    for i in 1:len
        xpos = index_range_zyx[3][i] - xcenter
        res[param_idx, i] = (xpos * xscale)^order_zyx[3] * coeff
    end
    return res
end

const term_0 = generate_term(center_index_zyx, index_range_zyx, solution.stride,
                             (0, 0, 0))
const term_0_flat = vec(term_0)

const coeff = Matrix{Float64}(undef, length(term_0_flat), length(ele_data) + 1)
coeff[:, 1] .= term_0_flat
for i in 1:length(ele_data)
    coeff[:, 1 + i] .= vec(ele_data[i])
end

function fit_term(model::Model, target, coeff, maxv, yz_weight)
    nvars = size(coeff, 2)
    @variable(model, -maxv .<= vars[1:nvars] .<= maxv)
    # No limit on the DC term
    delete_lower_bound(vars[1])
    delete_upper_bound(vars[1])
    err = coeff * vars .- vec(target)
    nvals = length(err)
    @variable(model, abserr[1:nvals])
    @constraint(model, .-abserr .<= err)
    @constraint(model, abserr .>= err)
    abserr = reshape(abserr, size(target))

    @variable(model, maxerr)

    npoints = size(target, 2)
    point_center = (1 + npoints) / 2
    max_offset = point_center - 1
    C = log(2) / max_offset^2
    for i in 1:size(target, 2)
        x_weight = exp(-(i - point_center)^2 * C)
        for j in 1:6
            e = abserr[j, i] * x_weight * yz_weight[j]
            @constraint(model, -maxerr <= e)
            @constraint(model, maxerr >= e)
        end
    end
    @objective(model, Min, maxerr)
    optimize!(model)
    @show model
    return @show value.(vars)
end

const term_x1 = generate_term(center_index_zyx, index_range_zyx, solution.stride,
                              (0, 0, 1))
const term_x2_raw = generate_term(center_index_zyx, index_range_zyx, solution.stride,
                                  (0, 0, 2))
const term_y2_raw = generate_term(center_index_zyx, index_range_zyx, solution.stride,
                                  (0, 2, 0))
const term_z2_raw = generate_term(center_index_zyx, index_range_zyx, solution.stride,
                                  (2, 0, 0))

const term_x2 = term_x2_raw .- (term_y2_raw .+ term_z2_raw) ./ 2
const term_zz = term_y2_raw .- term_z2_raw

const term_x1_flat = vec(term_x1)
const term_x2_flat = vec(term_x2)
const term_zz_flat = vec(term_zz)

const voltages_x1 = coeff \ term_x1_flat
const voltages_x2 = coeff \ term_x2_flat
const voltages_zz = coeff \ term_zz_flat

const voltages_x1_2 = fit_term(Model(Ipopt.Optimizer), term_x1, coeff, 0.3,
                               (4, 3, 3, 1, 3, 1))
const voltages_x2_2 = fit_term(Model(Ipopt.Optimizer), term_x2, coeff, 10.0,
                               (40, 3, 3, 2000, 300, 2000))
const voltages_zz_2 = fit_term(Model(Ipopt.Optimizer), term_zz, coeff, 10.0,
                               (40, 3, 3, 2000, 300, 2000))

function get_nth_part(data, idx)
    return data[idx:6:end]
end

const term_names = ["1", "y", "z", "y2", "yz", "z2"]

@show extrema(voltages_x1[2:end])
true_x1 = coeff * voltages_x1
err_x1 = term_x1_flat .- true_x1
@show extrema(err_x1)

@show extrema(voltages_x1_2[2:end])
true_x1_2 = coeff * voltages_x1_2
err_x1_2 = term_x1_flat .- true_x1_2
@show extrema(err_x1_2)

@show extrema(voltages_x2[2:end])
true_x2 = coeff * voltages_x2
err_x2 = term_x2_flat .- true_x2
@show extrema(err_x2)

@show extrema(voltages_x2_2[2:end])
true_x2_2 = coeff * voltages_x2_2
err_x2_2 = term_x2_flat .- true_x2_2
@show extrema(err_x2_2)

@show extrema(voltages_zz[2:end])
true_zz = coeff * voltages_zz
err_zz = term_zz_flat .- true_zz
@show extrema(err_zz)

@show extrema(voltages_zz_2[2:end])
true_zz_2 = coeff * voltages_zz_2
err_zz_2 = term_zz_flat .- true_zz_2
@show extrema(err_zz_2)


# for i in 1:6
#     figure()
#     plot(get_nth_part(term_x1_flat, i), label="X1 tgt")
#     plot(get_nth_part(true_x1, i), label="X1")
#     plot(get_nth_part(true_x1_2, i), label="X1_2")
#     legend()
#     title(term_names[i])
# end

# for i in 1:6
#     figure()
#     plot(get_nth_part(term_x2_flat, i), label="X2 tgt")
#     plot(get_nth_part(true_x2, i), label="X2")
#     plot(get_nth_part(true_x2_2, i), label="X2_2")
#     legend()
#     title(term_names[i])
# end

for i in 1:6
    figure()
    plot(get_nth_part(term_zz_flat, i), label="ZZ tgt")
    plot(get_nth_part(true_zz, i), label="ZZ")
    plot(get_nth_part(true_zz_2, i), label="ZZ_2")
    legend()
    title(term_names[i])
end

# for i in 1:6
#     figure()
#     plot(get_nth_part(err_x1, i), label="X1")
#     plot(get_nth_part(err_x1_2, i), label="X1_2")
#     plot(get_nth_part(err_x2, i), label="X2")
#     plot(get_nth_part(err_x2_2, i), label="X2_2")
#     legend()
#     title(term_names[i])
# end

NaCsPlot.maybe_show()
