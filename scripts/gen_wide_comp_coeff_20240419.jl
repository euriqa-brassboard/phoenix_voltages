#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
import PhoenixVoltages.Fitting
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using MAT
using LinearAlgebra

# Loading input data
const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202310.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const center_fit_cache = Solutions.compensate_fitter3(solution)

# Output prefix
const prefix = joinpath(@__DIR__, "../data/wide_comp_coeff_20240419")

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
    eles::Vector{Int}
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

    ncenter_coeff = length(ele_data[1].center_data)
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
    return AllTermsData(eles, ele_data, center_x0,
                        center_B, slice_coeff, slice_terms)
end

function get_wide_comp_coeff(xpos_um)
    @show xpos_um

    all_term_data = AllTermsData(center_fit_cache, xpos_um)

    return Dict("electrodes"=>all_term_data.eles,
                "solution"=>all_term_data.center_x0,
                "free_solution"=>all_term_data.center_B,
                "slice_coeff"=>all_term_data.slice_coeff,
                "slice_terms"=>all_term_data.slice_terms,
                "ele_datas"=>[Dict("slice_data"=>ele_data.slice_data,
                                   "center_data"=>ele_data.center_data)
                              for ele_data in all_term_data.ele_data],
                "xpos_um"=>xpos_um)
end

const xpos_ums = -350:350
matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", solution.electrode_names)
end
mkpath("$(prefix)")
for xpos_um in xpos_ums
    matopen("$(prefix)/$(xpos_um).mat", "w") do mat
        write(mat, "data", @time(get_wide_comp_coeff(xpos_um)))
    end
end
