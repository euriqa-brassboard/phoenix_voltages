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
import Ipopt
using JuMP

const centers = Solutions.CenterTracker()
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)

struct Frame
    target::Vector{Float64}
    weight::Vector{Float64}
    freedom::Vector{Tuple{Vector{Float64},NTuple{2,Float64}}}
    electrodes::Vector{Int}
    electrode_potentials::Matrix{Float64}
end

struct TrapModel
    model::Model
    # Parameters to optimize (electrode voltages)
    vs::Vector{Vector{VariableRef}}

    # Cost functions for each frame
    # Maximum diviation from target potential
    maxdiffs::Vector{VariableRef}
    # Maximum voltage
    maxvs::Vector{VariableRef}
    # Maximum voltage difference with last frame
    maxvdiffs::Vector{VariableRef}

    # Global cost functions
    maxmaxv::VariableRef
    maxmaxvdiff::VariableRef

    function TrapModel()
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "max_iter", 30000)
        set_optimizer_attribute(model, "print_level", 5)
        return new(model, Vector{VariableRef}[],
                   VariableRef[], VariableRef[], VariableRef[],
                   @variable(model), @variable(model))
    end
end

struct TrapWeights
    maxdiffs::Float64
    maxvs::Float64
    maxvdiffs::Float64
    maxmaxv::Float64
    maxmaxvdiff::Float64
end

mutable struct TrapModelBuilder
    model::TrapModel
    nframes::Int
    prev_vmap::Dict{Int,VariableRef}
    TrapModelBuilder() = new(TrapModel(), 0)
end

function add_frame!(builder::TrapModelBuilder, frame::Frame)
    tm = builder.model
    model = tm.model
    neles = length(frame.electrodes)
    @assert size(frame.electrode_potentials, 2) == neles

    target = frame.target
    for flex in frame.freedom
        lb, ub = flex[2]
        if (lb == 0 && ub == 0) || (lb > ub)
            continue
        end
        c = @variable(model)
        target = @expression(model, target .+ c .* flex[1])
        if isfinite(lb)
            @constraint(model, c >= lb)
        end
        if isfinite(ub)
            @constraint(model, c <= ub)
        end
    end

    v = @variable(model, [1:neles])
    push!(tm.vs, v)

    all_diff = @expression(model, (frame.electrode_potentials * v .- target) .* frame.weight)
    maxdiff = @variable(model)
    @constraint(model, maxdiff .>= all_diff)
    @constraint(model, maxdiff .>= .-all_diff)
    push!(tm.maxdiffs, maxdiff)

    maxv = @variable(model)
    @constraint(model, maxv .>= v)
    @constraint(model, maxv .>= .-v)
    push!(tm.maxvs, maxv)

    vmap = Dict(zip(frame.electrodes, v))
    if isdefined(builder, :prev_vmap)
        maxvdiff = @variable(model)
        prev_vmap = builder.prev_vmap
        for k in union(keys(vmap), keys(prev_vmap))
            v1 = get(vmap, k, 0)
            v2 = get(prev_vmap, k, 0)
            @constraint(model, maxvdiff >= v1 - v2)
            @constraint(model, maxvdiff >= v2 - v1)
        end
        push!(tm.maxvdiffs, maxvdiff)
        @constraint(model, tm.maxmaxvdiff >= maxvdiff)
    end
    builder.prev_vmap = vmap

    @constraint(model, tm.maxmaxv >= maxv)
    builder.nframes += 1

    return builder
end

function finalize_trap_model!(builder::TrapModelBuilder, weights::TrapWeights)
    tm = builder.model
    model = tm.model

    obj = sum(tm.maxdiffs) * weights.maxdiffs +
        sum(tm.maxvs) * weights.maxvs
    if weights.maxvdiffs > 0
        obj += sum(tm.maxvdiffs) * weights.maxvdiffs
    end
    if weights.maxmaxv > 0
        obj += sum(tm.maxmaxv) * weights.maxmaxv * builder.nframes
    end
    if weights.maxmaxvdiff > 0
        obj += sum(tm.maxmaxvdiff) * weights.maxmaxvdiff * builder.nframes
    end
    @objective(model, Min, obj)
    return builder
end

function optimize_trap_model!(builder::TrapModelBuilder)
    JuMP.optimize!(builder.model.model)
    return [[value(v) for v in v] for v in builder.model.vs]
end

struct MultiFitterCache
    potential::Potential
    fit_caches::Dict{Int,Potentials.FitCache}
    function MultiFitterCache(potential::Potential)
        return new(potential, Dict{Int,Potentials.FitCache}())
    end
end

Base.get(cache::MultiFitterCache, xsize) = get!(cache.fit_caches, xsize) do
    if xsize > 50
        order = 8
    elseif xsize > 25
        order = 6
    else
        order = 4
    end
    fitter = Fitting.PolyFitter(2, 2, order, sizes=(5, 5, xsize))
    return Potentials.FitCache(fitter, cache.potential)
end
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
const center_x2 = 1.0

const move_relax_start_um = -250
const move_relax_end_um = -100
function get_move_x2(xpos_um)
    r = get_ratio(xpos_um, move_relax_start_um, move_relax_end_um)
    return interpolate(r, move_x2_init, center_x2)
end

const stride_um = solution.stride .* 1000

# Size of the longer side
# We are kind of directly assuming that the grid size is 1um here.
# Too lazy to apply the correction factor
@assert all(stride_um .== 1)

const move_fit_init_size = 14
const move_fit_min_size = 6
const center_fit_init_size = 31
const center_fit_min_size = 14

# Left/right assymetry
const max_fit_assym_ratio = 2.5
# Gap / (move total fit range + center total fit range)
const min_fit_gap_ratio = 1.5

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
const move_fit_final_mid_size = 0
const move_fit_final_total_size = move_fit_final_out_size + move_fit_final_mid_size
const center_fit_final_out_size = center_fit_init_size
const center_fit_final_mid_size = 0
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
        r = get_ratio(xpos_um, fit_shrink_end_um, 0)
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
    return round(Int, center_pos_idx[1] - move_fit_final_out_size):round(Int, center_pos_idx[1] + move_fit_final_out_size)
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

@enum(TermIndex, C=1, DX=2, DY=3, DZ=4, XY=5, YZ=6, ZX=7, Z2=8, X2=9)

function term_block(term, scale)
    block = zeros(Int(typemax(TermIndex)))
    block[Int(term)] = scale
    return block
end

c_block(scale=1) = term_block(C, scale)
dx_block(scale=1) = term_block(DX, scale)
dy_block(scale=1) = term_block(DY, scale)
dz_block(scale=1) = term_block(DZ, scale)

xy_block(scale=1) = term_block(XY, scale)
yz_block(scale=1) = term_block(YZ, scale)
zx_block(scale=1) = term_block(ZX, scale)

z2_block(scale=1) = term_block(Z2, scale)
x2_block(scale=1) = term_block(X2, scale)

function weight_block()
    block = zeros(Int(typemax(TermIndex)))
    block[Int(C)] = 1
    block[Int(DX)] = 1
    block[Int(DY)] = 1
    block[Int(DZ)] = 1

    block[Int(XY)] = 0.7
    block[Int(YZ)] = 0.7
    block[Int(ZX)] = 0.7

    block[Int(Z2)] = 0.7
    block[Int(X2)] = 1

    return block
end

function potential_block(ele, center, xrange)
    xleft = first(xrange)
    xright = last(xrange)
    xsize = xright - xleft + 1
    fit_cache = get(multi_fit_cache, xsize)
    center_r = (center[3], center[2], center[1])
    fit_center_r = (center[3], center[2], (xright + xleft) / 2)

    fit = get(fit_cache, ele, center_r, fit_center=fit_center_r)
    terms = Solutions.get_compensate_terms1(fit, stride_um)

    block = zeros(Int(typemax(TermIndex)))

    block[Int(C)] = fit[0, 0, 0]
    block[Int(DX)] = terms.dx
    block[Int(DY)] = terms.dy
    block[Int(DZ)] = terms.dz

    block[Int(XY)] = terms.xy
    block[Int(YZ)] = terms.yz
    block[Int(ZX)] = terms.zx

    block[Int(Z2)] = terms.z2
    block[Int(X2)] = terms.x2

    return block
end

const max_shape_relax_factor = 1
const shape_relax_start_um = -2900
const shape_relax_max_um = -100
const shape_relax_end_um = center_um
function get_shape_relax_factor(xpos_um)
    if xpos_um <= shape_relax_max_um
        r = get_ratio(xpos_um, shape_relax_start_um, shape_relax_max_um)
        return interpolate(r, 0, max_shape_relax_factor)
    else
        r = get_ratio(xpos_um, shape_relax_max_um, shape_relax_end_um)
        return interpolate(r, max_shape_relax_factor, 0)
    end
end

const dc_level_limit_max = 2.5
const dc_level_limit_start_um = -250
const dc_level_limit_end_um = 0
function get_dc_level_limit(xpos_um)
    if xpos_um <= dc_level_limit_start_um
        return Inf
    end
    r = get_ratio(xpos_um, dc_level_limit_start_um, dc_level_limit_end_um)
    return interpolate(r, dc_level_limit_max, 0)
end

function construct_frame(eles, target, weight, freedom, ranges, centers)
    neles = length(eles)
    _target = flatten_blocks(target)
    npoints = length(_target)
    _weight = flatten_blocks(weight)
    @assert length(_weight) == npoints
    _freedom = [(flatten_blocks(f[1]), (f[2], f[3])) for f in freedom]
    electrode_potentials = Matrix{Float64}(undef, npoints, neles)
    for i in 1:neles
        ele = eles[i]
        electrode_potentials[:, i] = flatten_blocks(potential_block.(ele, centers,
                                                                     ranges))
    end
    return Frame(_target, _weight, _freedom, eles, electrode_potentials)
end

function create_frame0(eles)
    center_xrange = get_xranges0()

    target_x2 = x2_block(center_x2)
    target_yz = xy_block(center_x2 / 2)
    target = target_x2 .+ target_yz

    weight = weight_block()

    freedom = Tuple{Vector{Float64},Float64,Float64}[]
    push!(freedom, (c_block(), -Inf, Inf))

    return construct_frame(eles, target, weight, freedom,
                           (center_xrange,), (center_pos_idx,))
end

function create_frame1(eles, xpos_um)
    move_pos_idx = get_rf_center(xpos_um)
    move_xrange, center_xrange = get_xranges1(xpos_um, move_pos_idx)

    shape_relax = get_shape_relax_factor(xpos_um)

    center_yz = center_x2 / 2
    center_target_x2 = x2_block(center_x2)
    center_target_yz = xy_block(center_yz)
    center_target = center_target_x2 .+ center_target_yz

    move_x2 = get_move_x2(xpos_um)
    move_yz = move_x2 / 2
    move_target_x2 = x2_block(move_x2)
    move_target_yz = xy_block(move_yz)
    move_target = move_target_x2 .+ move_target_yz
    target = center_target, move_target

    center_weight = weight_block()
    move_weight = weight_block()
    weight = center_weight, move_weight

    center_zero = c_block(0)
    move_zero = c_block(0)
    center_one = c_block(1)
    move_one = c_block(1)

    freedom = Tuple{NTuple{2,Vector{Float64}},Float64,Float64}[]
    push!(freedom, ((xy_block(), move_zero), -Inf, Inf))
    push!(freedom, ((yz_block(), move_zero),
                    -center_yz * shape_relax / 3, center_yz * shape_relax / 3))
    push!(freedom, ((z2_block(), move_zero),
                    -center_yz * shape_relax / 3, center_yz * shape_relax / 3))
    push!(freedom, ((x2_block(), move_zero),
                    -center_x2 * shape_relax / 5, center_x2 * shape_relax / 5))
    push!(freedom, ((zx_block(), move_zero),
                    -center_x2 * shape_relax / 5, center_x2 * shape_relax / 5))

    push!(freedom, ((center_zero, xy_block()), -Inf, Inf))
    push!(freedom, ((center_zero, yz_block()),
                    -move_yz * shape_relax / 3, move_yz * shape_relax / 3))
    push!(freedom, ((center_zero, z2_block()),
                    -move_yz * shape_relax / 3, move_yz * shape_relax / 3))
    push!(freedom, ((center_zero, x2_block()),
                    -move_x2 * shape_relax / 4, move_x2 * shape_relax / 4))
    push!(freedom, ((center_zero, zx_block()), -Inf, Inf))

    push!(freedom, ((center_one, move_one), -Inf, Inf))

    dc_level_limit = get_dc_level_limit(xpos_um)
    push!(freedom, ((center_one, .-move_one), -dc_level_limit, dc_level_limit))

    return construct_frame(eles, target, weight, freedom,
                           (center_xrange, move_xrange),
                           (center_pos_idx, move_pos_idx))
end

function create_frame(xpos_um)
    @show xpos_um
    eles = get_electrodes(xpos_um)
    if xpos_um == 0
        return eles, create_frame0(eles)
    else
        return eles, create_frame1(eles, xpos_um)
    end
end

const prefix = joinpath(@__DIR__, "../data/merge_data_20221001")

const xpos_ums = loading_um:center_um
const builder = TrapModelBuilder()

const results = Dict{String,Any}[]
for xpos_um in xpos_ums
    @time begin
        eles, frame = create_frame(xpos_um)
        add_frame!(builder, frame)
        push!(results, Dict("electrodes"=>eles, "xpos_um"=>xpos_um))
    end
end
@time finalize_trap_model!(builder, TrapWeights(1, 0.1, 0.1, 1, 1))
for (rd, vals) in zip(results, optimize_trap_model!(builder))
    rd["voltages"] = vals
end

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", solution.electrode_names)
    write(mat, "transfer_solutions", results)
end
