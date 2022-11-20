#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
using PhoenixVoltages.Solutions
using PhoenixVoltages.Solutions: l_unit, V_unit, l_unit_um, V_unit_uV
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

struct Frame
    solution::Vector{Float64}
    limited_solution::Vector{Vector{Float64}}
    limits::Vector{NTuple{2,Float64}}
    electrodes::Vector{Int}
    constraints::Vector{Tuple{Vector{Float64},Float64,Float64}}
    function Frame(target::Vector, freedom::Vector, electrodes::Vector{Int},
                   electrode_potentials::Matrix)
        neles = length(electrodes)
        @assert size(electrode_potentials, 2) == neles

        limited_solution = Vector{Float64}[]
        limits = NTuple{2,Float64}[]
        for flex in freedom
            lb, ub = flex[2]
            if lb > ub
                continue
            end
            if lb == ub
                if lb != 0
                    target = target .+ lb .* flex[1]
                end
                continue
            end
            push!(limited_solution, electrode_potentials \ flex[1])
            push!(limits, flex[2])
        end
        v = electrode_potentials \ target
        coeff = electrode_potentials
        B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]
        nv, nt = size(B)
        @assert length(v) == nv
        for i in 1:nt
            push!(limited_solution, B[:, i])
            push!(limits, (-Inf, Inf))
        end

        return new(v, limited_solution, limits, electrodes,
                   Tuple{Vector{Float64},Float64,Float64}[])
    end
end

struct MultiFitterCache
    potential::Potential
    fit_caches::Dict{Int,Potentials.FitCache}
    function MultiFitterCache(potential::Potential)
        return new(potential, Dict{Int,Potentials.FitCache}())
    end
end

function Base.get(cache::MultiFitterCache, xsize)
    if xsize in keys(cache.fit_caches)
        return cache.fit_caches[xsize]
    end
    if length(cache.fit_caches) > 10
        empty!(cache.fit_caches)
    end
    if xsize > 150
        order = 16
    elseif xsize > 125
        order = 14
    elseif xsize > 100
        order = 12
    elseif xsize > 75
        order = 10
    elseif xsize > 50
        order = 8
    elseif xsize > 25
        order = 6
    else
        order = 4
    end
    fitter = Fitting.PolyFitter(2, 2, order, sizes=(5, 5, xsize))
    res = Potentials.FitCache(fitter, cache.potential)
    cache.fit_caches[xsize] = res
    return res
end
Base.empty!(cache::MultiFitterCache) = empty!(cache.fit_caches)
const multi_fit_cache = MultiFitterCache(solution)

flatten_blocks(blocks::Union{NTuple{N,A} where N,Vector{A}} where A <: Array{Float64}) =
    reduce(vcat, vec.(blocks))
flatten_blocks(block::Array{Float64}) = vec(block)

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_ratio(x, x0, x1)::Float64
    if x0 < x1
        if x <= x0
            return 0
        elseif x >= x1
            return 1
        end
    else
        if x <= x1
            return 1
        elseif x >= x0
            return 0
        end
    end
    return (x - x0) / (x1 - x0)
end
function interpolate(x, v0, v1)
    if x <= 0
        return v0
    elseif x >= 1
        return v1
    end
    return x * v1 + (1 - x) * v0
end

const loading_um = -3045
const center_um = 0

const center_electrodes = Mappings.find_electrodes(solution.electrode_index,
                                                   center_um, min_num=20,
                                                   min_dist=350)
function get_electrodes(xpos_um)
    if abs(xpos_um - center_um) <= 100
        eles = center_electrodes
    else
        eles = Mappings.find_electrodes(solution.electrode_index,
                                        xpos_um, min_num=15, min_dist=250)
        union!(eles, center_electrodes)
    end
    return sort!(collect(eles))
end

const move_x2_init = 3.0
const move_x2_final = 1.0
const center_x2 = 1.0

const move_relax_start_um = -250
const move_relax_end_um = -160
function get_move_x2(xpos_um)
    r = get_ratio(xpos_um, move_relax_start_um, move_relax_end_um)
    return interpolate(r, move_x2_init, move_x2_final)
end

const stride_um = solution.stride .* 1000

# Size of the longer side
# We are kind of directly assuming that the grid size is 1um here.
# Too lazy to apply the correction factor
@assert all(stride_um .== 1)

const move_fit_init_size = 25
const move_fit_min_size = 25
const center_fit_init_size = 51
const center_fit_min_size = 27

# Left/right assymetry
const max_fit_assym_ratio = 2.5
# Gap / (move total fit range + center total fit range)
const min_fit_gap_ratio = 1.2

const move_fit_init_out_size = move_fit_init_size
const move_fit_init_mid_size = move_fit_init_size
const move_fit_init_total_size = move_fit_init_size * 2
const center_fit_init_out_size = center_fit_init_size
const center_fit_init_mid_size = center_fit_init_size
const center_fit_init_total_size = center_fit_init_size * 2

const move_fit_shrunk_out_size = move_fit_min_size
const move_fit_shrunk_mid_size = move_fit_min_size / max_fit_assym_ratio
const move_fit_shrunk_total_size = move_fit_shrunk_out_size + move_fit_shrunk_mid_size
const center_fit_shrunk_out_size = center_fit_min_size
const center_fit_shrunk_mid_size = center_fit_min_size / max_fit_assym_ratio
const center_fit_shrunk_total_size =
    center_fit_shrunk_out_size + center_fit_shrunk_mid_size

const move_fit_final_out_size = center_fit_init_size
const move_fit_final_mid_size = -3
const move_fit_final_total_size = move_fit_final_out_size + move_fit_final_mid_size
const center_fit_final_out_size = center_fit_init_size
const center_fit_final_mid_size = -3
const center_fit_final_total_size =
    center_fit_final_out_size + center_fit_final_mid_size

# When the initial fitting size first start to touch each other
const fit_shrink_start_um = -((move_fit_init_mid_size + center_fit_init_mid_size) +
    (move_fit_init_total_size + center_fit_init_total_size) * min_fit_gap_ratio)
# When we reach the limit of the shrinking
const fit_shrink_end_um = -((move_fit_shrunk_mid_size + center_fit_shrunk_mid_size) +
    (move_fit_shrunk_total_size + center_fit_shrunk_total_size) * min_fit_gap_ratio)

function get_xspans(xpos_um)
    # The caller should handle the center case explicitly
    @assert xpos_um <= -1
    if xpos_um <= fit_shrink_end_um
        r = get_ratio(xpos_um, fit_shrink_start_um, fit_shrink_end_um)
        move_out = interpolate(r, move_fit_init_out_size, move_fit_shrunk_out_size)
        move_mid = interpolate(r, move_fit_init_mid_size, move_fit_shrunk_mid_size)
        center_out = interpolate(r, center_fit_init_out_size,
                                 center_fit_shrunk_out_size)
        center_mid = interpolate(r, center_fit_init_mid_size,
                                 center_fit_shrunk_mid_size)
    else
        r = get_ratio(xpos_um, fit_shrink_end_um, center_um)
        move_out = interpolate(r, move_fit_shrunk_out_size, move_fit_final_out_size)
        move_mid = interpolate(r, move_fit_shrunk_mid_size, move_fit_final_mid_size)
        center_out = interpolate(r, center_fit_shrunk_out_size,
                                 center_fit_final_out_size)
        center_mid = interpolate(r, center_fit_shrunk_mid_size,
                                 center_fit_final_mid_size)
    end
    move_left = -move_out
    move_right = move_mid
    center_left = -center_mid
    center_right = center_out

    return (move_left, move_right), (center_left, center_right)
end

const center_pos_idx = get_rf_center(center_um)
function get_xranges0()
    return round(Int, center_pos_idx[1] - 136):round(Int, center_pos_idx[1] + 136)
end
function get_xranges1(xpos_um, move_pos_idx)
    (move_left, move_right), (center_left, center_right) = get_xspans(xpos_um)

    move_left = round(Int, move_pos_idx[1] + move_left)
    move_right = round(Int, move_pos_idx[1] + move_right)

    center_left = round(Int, center_pos_idx[1] + center_left)
    center_right = round(Int, center_pos_idx[1] + center_right)
    @assert center_left > move_right
    return (move_left:move_right, center_left:center_right)
end

@enum(TermIndex, C=1, DX=2, DY=3, DZ=4, XY=5, YZ=6, ZX=7, Z2=8, X2=9, X3=10, X4=11,
      X5=12, X6=13, X7=14, X8=15)
const NTerms1 = 11
const NTerms2 = 15

function zero_block(version)
    @assert version == 1 || version == 2
    return zeros(version == 1 ? Int(NTerms1) : Int(NTerms2))
end

function term_block(term, scale, version)
    block = zero_block(version)
    block[Int(term)] = scale
    return block
end

c_block(scale=1, version=1) = term_block(C, scale, version)
dx_block(scale=1, version=1) = term_block(DX, scale, version)
dy_block(scale=1, version=1) = term_block(DY, scale, version)
dz_block(scale=1, version=1) = term_block(DZ, scale, version)

xy_block(scale=1, version=1) = term_block(XY, scale, version)
yz_block(scale=1, version=1) = term_block(YZ, scale, version)
zx_block(scale=1, version=1) = term_block(ZX, scale, version)

z2_block(scale=1, version=1) = term_block(Z2, scale, version)
x2_block(scale=1, version=1) = term_block(X2, scale, version)
x3_block(scale=1, version=1) = term_block(X3, scale, version)
x4_block(scale=1, version=1) = term_block(X4, scale, version)

x5_block(scale=1, version=1) = term_block(X5, scale, version)
x6_block(scale=1, version=1) = term_block(X6, scale, version)
x7_block(scale=1, version=1) = term_block(X7, scale, version)
x8_block(scale=1, version=1) = term_block(X8, scale, version)

function potential_block(ele, center, xrange, version=1)
    xleft = first(xrange)
    xright = last(xrange)
    xsize = xright - xleft + 1
    fit_cache = get(multi_fit_cache, xsize)
    center_r = (center[3], center[2], center[1])
    fit_center_r = (center[3], center[2], (xright + xleft) / 2)

    fit = get(fit_cache, ele, center_r, fit_center=fit_center_r)
    terms = Solutions.get_compensate_terms1(fit, stride_um)

    block = zero_block(version)

    block[Int(C)] = fit[0, 0, 0]
    block[Int(DX)] = terms.dx
    block[Int(DY)] = terms.dy
    block[Int(DZ)] = terms.dz

    block[Int(XY)] = terms.xy
    block[Int(YZ)] = terms.yz
    block[Int(ZX)] = terms.zx

    block[Int(Z2)] = terms.z2
    block[Int(X2)] = terms.x2
    block[Int(X3)] = terms.x3
    block[Int(X4)] = terms.x4

    if version >= 2
        block[Int(X5)] = fit[0, 0, 5]
        block[Int(X6)] = fit[0, 0, 6]
        block[Int(X7)] = fit[0, 0, 7]
        block[Int(X8)] = fit[0, 0, 8]
    end

    return block
end

const max_shape_relax_factor = 1
const shape_relax_start_um = -2900
const shape_relax_max_um = -30
const shape_relax_end_um = center_um
function get_shape_relax_factor(xpos_um)
    if xpos_um <= shape_relax_max_um
        r = get_ratio(xpos_um, shape_relax_start_um, shape_relax_max_um)
        return interpolate(r, 0, max_shape_relax_factor)
    elseif xpos_um == 0
        return 0.0
    else
        r = get_ratio(xpos_um, shape_relax_max_um, shape_relax_end_um)
        r = 1 - sqrt(1 - r)
        return interpolate(r, max_shape_relax_factor, 0)
    end
end

const dc_level_limit_max = 0.03
const dc_level_limit_start_um = -500
const dc_level_limit_end_um = -240
function get_dc_level_limit(xpos_um)
    if xpos_um <= dc_level_limit_start_um
        return Inf
    end
    r = get_ratio(xpos_um, dc_level_limit_start_um, dc_level_limit_end_um)
    return interpolate(r, dc_level_limit_max, 0)
end

function get_barriers(xpos_um)
    barrier_dist = 175
    barrier_max_height = 1.0

    barriers = NTuple{3,Float64}[(xpos_um - barrier_dist, xpos_um, barrier_max_height),
                                 (center_um + barrier_dist, center_um, barrier_max_height)]

    r = get_ratio(center_um - xpos_um, barrier_dist * 3.5, barrier_dist * 7)
    mid_height = interpolate(r, 0.0, barrier_max_height)
    if mid_height > 0
        push!(barriers, (xpos_um + barrier_dist, xpos_um, mid_height))
        push!(barriers, (center_um - barrier_dist, center_um, mid_height))
    end
    return barriers
end

function construct_frame(eles, target, freedom, ranges, centers, version=1)
    neles = length(eles)
    _target = flatten_blocks(target)
    npoints = length(_target)
    _freedom = [(flatten_blocks(f[1]), (f[2], f[3])) for f in freedom]
    electrode_potentials = Matrix{Float64}(undef, npoints, neles)
    for i in 1:neles
        ele = eles[i]
        electrode_potentials[:, i] = flatten_blocks(potential_block.(ele, centers,
                                                                     ranges, version))
    end
    return Frame(_target, _freedom, eles, electrode_potentials)
end

const frame_cutoff_um = 0
# const frame_cutoff_um = -10

# function compute_start_o4_coeff0()
#     move_x2 = get_move_x2(frame_cutoff_um)
#     move_pos_um = frame_cutoff_um - center_um
#     α = move_x2 / center_x2

#     a = (α + 1) / move_pos_um^2 * center_x2 * 6
#     b = -(α + 2) / move_pos_um * center_x2 * 2
#     c = center_x2
#     # The functional form initially is (a/24 * x^4 + b/6 * x^3 + c/2 x^2)
#     # in EURIQA unit
#     barrier_um = (-(b / 2) + sqrt((b / 2)^2 - 4 * (a / 6) * c)) / (2 * (a / 6))
#     # We'll translate the coordinate to be centered around the barrier @ barrier_um < 0
#     a1 = a
#     b1 = @evalpoly(barrier_um, b, a)
#     c1 = @evalpoly(barrier_um, c, b, a / 2)
#     # d1 = @evalpoly(barrier_um, 0, c, b / 2, a / 6)
#     return a1, b1, c1, barrier_um
# end
# const start_o4_coeff0 = compute_start_o4_coeff0()
# const final_o4_coeff0 = (0.0, 0.0, center_x2, 0.0)

# function compute_o4_coeff0(xpos_um)
#     r = get_ratio(xpos_um, frame_cutoff_um, center_um)
#     a0, b0, c0, shift_um = interpolate.(r, start_o4_coeff0, final_o4_coeff0)
#     a = a0
#     b = @evalpoly(-shift_um, b0, a0)
#     c = @evalpoly(-shift_um, c0, b0, a0 / 2)
#     d = @evalpoly(-shift_um, 0, c0, b0 / 2, a0 / 6)
#     @show (a, b, c, d / (V_unit_uV / l_unit_um))
#     return (a, b, c, d / (V_unit_uV / l_unit_um))
# end

# function create_frame0(eles, xpos_um)
#     center_xrange = get_xranges0()

#     target = zero_block(2)
#     target[Int(X4)], target[Int(X3)], target[Int(X2)], target[Int(DX)] =
#         compute_o4_coeff0(xpos_um)

#     shape_relax = get_shape_relax_factor(xpos_um)

#     freedom = Tuple{Vector{Float64},Float64,Float64}[]
#     push!(freedom, (dx_block(1, 2), -3000 * shape_relax, 3000 * shape_relax))
#     push!(freedom, (dy_block(1, 2), -10 * shape_relax, 10 * shape_relax))
#     push!(freedom, (dz_block(1, 2), -10 * shape_relax, 10 * shape_relax))
#     push!(freedom, (xy_block(1, 2),
#                     -center_x2 * shape_relax / 5, center_x2 * shape_relax / 5))
#     push!(freedom, (yz_block(1, 2),
#                     -center_x2 * shape_relax / 8, center_x2 * shape_relax / 8))
#     push!(freedom, (zx_block(1, 2),
#                     -center_x2 * shape_relax / 5, center_x2 * shape_relax / 5))
#     push!(freedom, (z2_block(1, 2),
#                     -center_x2 * shape_relax / 8, center_x2 * shape_relax / 8))
#     push!(freedom, (x2_block(1, 2),
#                     -center_x2 * shape_relax / 5, center_x2 * shape_relax / 5))
#     push!(freedom, (c_block(1, 2), -Inf, Inf))

#     return construct_frame(eles, target, freedom,
#                            (center_xrange,), (center_pos_idx,), 2)
# end

function create_frame0(eles, xpos_um)
    @assert xpos_um == 0
    center_xrange = get_xranges0()

    center_yz = center_x2 / 2
    target_x2 = x2_block(center_x2, 2)
    target_yz = yz_block(center_yz, 2)
    target = target_x2 .+ target_yz

    freedom = Tuple{Vector{Float64},Float64,Float64}[]
    push!(freedom, (c_block(1, 2), -Inf, Inf))
    push!(freedom, (x6_block(1, 2), 0, Inf))
    push!(freedom, (x8_block(1, 2), 0, Inf))
    return construct_frame(eles, target, freedom,
                           (center_xrange,), (center_pos_idx,), 2)
end

function create_frame1(eles, xpos_um)
    move_pos_idx = get_rf_center(xpos_um)
    move_xrange, center_xrange = get_xranges1(xpos_um, move_pos_idx)

    shape_relax = get_shape_relax_factor(xpos_um)

    center_yz = center_x2 / 2
    center_target_x2 = x2_block(center_x2)
    center_target_yz = yz_block(center_yz)
    center_target = center_target_x2 .+ center_target_yz

    move_x2 = get_move_x2(xpos_um)
    move_yz = move_x2 / 2
    move_target_x2 = x2_block(move_x2)
    move_target_yz = yz_block(move_yz)
    move_target = move_target_x2 .+ move_target_yz
    target = center_target, move_target

    center_zero = c_block(0)
    move_zero = c_block(0)
    center_one = c_block(1)
    move_one = c_block(1)

    freedom = Tuple{NTuple{2,Vector{Float64}},Float64,Float64}[]
    push!(freedom, ((dx_block(), move_zero), -1250 * shape_relax, 1250 * shape_relax))
    push!(freedom, ((yz_block(), move_zero),
                    -center_yz * shape_relax / 4, center_yz * shape_relax / 4))
    push!(freedom, ((z2_block(), move_zero),
                    -center_yz * shape_relax / 4, center_yz * shape_relax / 4))
    push!(freedom, ((x2_block(), move_zero),
                    -center_x2 * shape_relax / 5, center_x2 * shape_relax / 5))
    push!(freedom, ((x3_block(), move_zero),
                    -center_x2 * shape_relax * 0.2, center_x2 * shape_relax * 1))
    push!(freedom, ((x4_block(), move_zero), 0, Inf))

    push!(freedom, ((center_zero, dx_block()), -1250 * shape_relax, 1250 * shape_relax))
    push!(freedom, ((center_zero, yz_block()),
                    -move_yz * shape_relax / 3, move_yz * shape_relax / 3))
    push!(freedom, ((center_zero, zx_block()),
                    -move_yz * shape_relax / 3, move_yz * shape_relax / 3))
    push!(freedom, ((center_zero, z2_block()),
                    -move_yz * shape_relax / 3, move_yz * shape_relax / 3))
    push!(freedom, ((center_zero, x2_block()),
                    -move_x2 * shape_relax / 8, move_x2 * shape_relax / 8))
    push!(freedom, ((center_zero, x3_block()),
                    -move_x2 * shape_relax * 1, move_x2 * shape_relax * 0.3))
    if xpos_um > -350
        push!(freedom, ((center_zero, x4_block()), 0, Inf))
    else
        push!(freedom, ((center_zero, x4_block()), -Inf, Inf))
    end

    push!(freedom, ((center_one, move_one), -Inf, Inf))

    dc_level_limit = get_dc_level_limit(xpos_um)
    if xpos_um == loading_um
        # Try balancing the DC level during loading
        dc_level_limit = 0.0
    end
    push!(freedom, ((center_one, .-move_one), -dc_level_limit, dc_level_limit))

    return construct_frame(eles, target, freedom,
                           (center_xrange, move_xrange),
                           (center_pos_idx, move_pos_idx))
end

function get_potential_diff_for_term(frame, potential, pos, base, term)
    pos_inds = get_rf_center(pos)
    pos_inds_r = (pos_inds[3], pos_inds[2], pos_inds[1])
    base_inds = get_rf_center(base)
    base_inds_r = (base_inds[3], base_inds[2], base_inds[1])
    fit_cache = get(multi_fit_cache, 9)
    res = 0.0
    for (e, v) in zip(frame.electrodes, term)
        res += v * (Potentials.get_single(fit_cache, e, pos_inds_r, (0, 0, 0)) -
            Potentials.get_single(fit_cache, e, base_inds_r, (0, 0, 0)))
    end
    return res
end

function add_barrier!(frame, xpos_um)
    barriers = get_barriers(xpos_um)
    terms = Iterators.flatten(((frame.solution,), frame.limited_solution))
    for (pos, base_pos, height) in barriers
        push!(frame.constraints, ([get_potential_diff_for_term(frame, solution, pos,
                                                               base_pos, term)
                                   for term in terms], height, Inf))
    end
end

function create_frame(xpos_um)
    @show xpos_um
    eles = get_electrodes(xpos_um)
    if xpos_um >= frame_cutoff_um
        frame = create_frame0(eles, xpos_um)
    else
        frame = create_frame1(eles, xpos_um)
    end
    add_barrier!(frame, xpos_um)
    return frame
end

const prefix = joinpath(@__DIR__, "../data/merge_coeff_20221118")

function dump_frame(xpos_um, frame)
    return Dict("electrodes"=>frame.electrodes,
                "solution"=>frame.solution,
                "limited_solution"=>frame.limited_solution,
                "limits"=>[v[i] for v in frame.limits, i in 1:2],
                "constraints_terms"=>[c[1] for c in frame.constraints],
                "constraints_limits"=>[c[i] for c in frame.constraints, i in 2:3],
                "xpos_um"=>xpos_um)
end

const xpos_ums = loading_um:center_um
# const xpos_ums = -500:center_um

const results = Dict{String,Any}[]
for xpos_um in xpos_ums
    @time begin
        frame = create_frame(xpos_um)
        push!(results, dump_frame(xpos_um, frame))
    end
end

matopen("$(prefix).mat", "w") do mat
    write(mat, "data", results)
    write(mat, "electrode_names", solution.electrode_names)
end
