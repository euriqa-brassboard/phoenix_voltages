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
    target::Vector{Float64}
    freedom::Vector{Tuple{Vector{Float64},NTuple{2,Float64}}}
    electrodes::Vector{Int}
    electrode_potentials::Matrix{Float64}
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

const loading_um = -3045.0
const center_um = 105.0
const parking_um = 455.0
const stride_um = solution.stride .* 1000
const xstride_um = stride_um[3]

const shuttling_steps = round(Int, abs(center_um - loading_um))
const merging_steps = round(Int, abs(parking_um - center_um))
const moveback_steps = merging_steps

mutable struct TrapSite
    xpos_um::Float64
    is_chain::Bool
    x2::Float64
    xidx_min::Int
    xidx_max::Int
    xorder::Int
    shape_relax::Float64
    level_min::Float64
    level_max::Float64
    electrodes::Set{Int}
    TrapSite(xpos_um, is_chain) = new(xpos_um, is_chain)
end

## Set position and type
# 0-based index
function sites_shuttling(idx)
    x = idx / shuttling_steps
    return [TrapSite(interpolate(x, loading_um, center_um), false),
            TrapSite(interpolate(x, center_um, parking_um), true)]
end
function sites_merging(idx)
    if idx == merging_steps
        return [TrapSite(center_um, true)]
    end
    x = idx / merging_steps
    return [TrapSite(center_um, false),
            TrapSite(interpolate(x, parking_um, center_um), true)]
end
function sites_moveback(idx)
    x = idx / moveback_steps
    return [TrapSite(interpolate(x, parking_um, center_um), true)]
end

function sites_gap(sites)
    if length(sites) == 1
        return -1
    end
    @assert length(sites) == 2
    return abs(sites[1].xpos_um - sites[2].xpos_um)
end

## Set trap shape
function set_x2!(sites)
    single_x2_init = 3.0
    single_x2_final = 1.0
    chain_x2 = 1.0

    relax_max_um = 250
    relax_min_um = 160

    dist = sites_gap(sites)

    if dist < 0
        @assert length(sites) == 1
        @assert sites[1].is_chain
        sites[1].x2 = chain_x2
        return
    end

    r = get_ratio(dist, relax_min_um, relax_max_um)
    single_x2 = interpolate(r, single_x2_final, single_x2_init)
    for site in sites
        site.x2 = site.is_chain ? chain_x2 : single_x2
    end
end

## Set fit range
function set_xrange!(site::TrapSite, left, right)
    xidx = Solutions.x_axis_to_index(solution, site.xpos_um / 1000)
    site.xidx_min = round(Int, xidx - left)
    site.xidx_max = round(Int, xidx + right)
    xsize = right + left
    # Intentionally smaller than the fitting one
    if xsize > 100
        site.xorder = 8
    elseif xsize > 50
        site.xorder = 6
    else
        site.xorder = 4
    end
    return
end

function set_xranges!(sites)
    single_out_size = 25
    chain_out_size = 51
    chain_max_size = 91
    inner_min_size = -3

    min_gap_ratio = 1.9
    min_gap_um = (single_out_size + chain_out_size) * min_gap_ratio * xstride_um

    dist = sites_gap(sites)
    if dist < 0
        @assert length(sites) == 1
        @assert sites[1].is_chain
        set_xrange!(sites[1], chain_max_size, chain_max_size)
        return
    end

    r = get_ratio(dist, 0, min_gap_um)
    chain_inner = interpolate(r, inner_min_size, chain_out_size)
    single_inner = interpolate(r, inner_min_size, single_out_size)
    single_outer = interpolate(r, single_out_size, chain_out_size)
    # Assume that single is on the left and is the first one in the vector...
    @assert !sites[1].is_chain
    @assert sites[2].is_chain
    @assert sites[1].xpos_um < sites[2].xpos_um
    for site in sites
        if site.is_chain
            set_xrange!(site, chain_inner, chain_out_size)
        else
            set_xrange!(site, single_outer, single_inner)
        end
    end
end

function set_shape_relax!(sites)
    max_shape_relax = 1.0
    relax_start_um = 800
    relax_max_um = 30
    relax_end_um = 0

    dist = sites_gap(sites)
    if dist < 0
        @assert length(sites) == 1
        @assert sites[1].is_chain
        sites[1].shape_relax = 0
        return
    end

    if dist >= relax_max_um
        r = get_ratio(dist, relax_max_um, relax_start_um)
        factor = interpolate(r, max_shape_relax, 0.0)
    else
        r = sqrt(get_ratio(dist, relax_end_um, relax_max_um))
        factor = interpolate(r, 0.0, max_shape_relax)
    end
    for site in sites
        site.shape_relax = factor
    end
end

function set_level_limit!(sites)
    level_limit_max = 0.08
    level_limit_start_um = 500
    level_limit_end_um = 240

    dist = sites_gap(sites)
    if dist < 0
        @assert length(sites) == 1
        @assert sites[1].is_chain
        sites[1].level_min = 0
        sites[1].level_max = 0
        return
    end

    if dist >= level_limit_start_um
        limit = Inf
    else
        r = get_ratio(dist, level_limit_end_um, level_limit_start_um)
        limit = interpolate(r, 0, level_limit_max)
    end

    for site in sites
        if site.is_chain
            site.level_min = 0
            site.level_max = 0
        else
            site.level_min = -limit
            site.level_max = limit
        end
    end
end

function set_electrodes!(sites)
    for site in sites
        if site.is_chain
            min_num = 20
            min_dist = 350
        else
            min_num = 15
            min_dist = 250
        end
        site.electrodes = Mappings.find_electrodes(solution.electrode_index,
                                                   site.xpos_um, min_num=min_num,
                                                   min_dist=min_dist)
    end
end

function get_electrodes(sites)
    eles = Set{Int}()
    for site in sites
        union!(eles, site.electrodes)
    end
    return sort!(collect(eles))
end

function populate_sites!(sites)
    set_x2!(sites)
    set_xranges!(sites)
    set_shape_relax!(sites)
    set_level_limit!(sites)
    set_electrodes!(sites)
end

@enum(BaseTermIndex, C=1, DX=2, DY=3, DZ=4, XY=5, YZ=6, ZX=7, Z2=8, X2=9)
const NBaseTerms = typemax(BaseTermIndex)

function zero_block(xorders)
    @assert xorders >= 2
    return zeros(Int(NBaseTerms) + xorders - 2)
end

function term_block(term, scale, xorders)
    block = zero_block(xorders)
    block[Int(term)] = scale
    return block
end

c_block(scale=1, xorders=4) = term_block(C, scale, xorders)
dx_block(scale=1, xorders=4) = term_block(DX, scale, xorders)
dy_block(scale=1, xorders=4) = term_block(DY, scale, xorders)
dz_block(scale=1, xorders=4) = term_block(DZ, scale, xorders)

xy_block(scale=1, xorders=4) = term_block(XY, scale, xorders)
yz_block(scale=1, xorders=4) = term_block(YZ, scale, xorders)
zx_block(scale=1, xorders=4) = term_block(ZX, scale, xorders)

z2_block(scale=1, xorders=4) = term_block(Z2, scale, xorders)
x2_block(scale=1, xorders=4) = term_block(X2, scale, xorders)
function xn_block(order, scale=1, xorders=4)
    return term_block(Int(X2) + order - 2, scale, xorders)
end

function potential_block(ele, center, xrange, xorders=4)
    xleft, xright = xrange
    xsize = xright - xleft + 1
    fit_cache = get(multi_fit_cache, xsize)
    center_r = (center[3], center[2], center[1])
    fit_center_r = (center[3], center[2], (xright + xleft) / 2)

    fit = get(fit_cache, ele, center_r, fit_center=fit_center_r)
    terms = Solutions.get_compensate_terms1(fit, stride_um)

    block = zero_block(xorders)

    block[Int(C)] = fit[0, 0, 0]
    block[Int(DX)] = terms.dx
    block[Int(DY)] = terms.dy
    block[Int(DZ)] = terms.dz

    block[Int(XY)] = terms.xy
    block[Int(YZ)] = terms.yz
    block[Int(ZX)] = terms.zx

    block[Int(Z2)] = terms.z2
    block[Int(X2)] = terms.x2
    for order in 3:xorders
        term = fit[0, 0, order]
        for i in 1:order
            term = term / xstride_um * l_unit_um * i
        end
        block[Int(X2) + order - 2] = term / V_unit
    end
    return block
end

function construct_frame(eles, target, freedom, ranges, centers, xorders)
    neles = length(eles)
    _target = flatten_blocks(target)
    npoints = length(_target)
    _freedom = [(flatten_blocks(f[1]), (f[2], f[3])) for f in freedom]
    electrode_potentials = Matrix{Float64}(undef, npoints, neles)
    for i in 1:neles
        ele = eles[i]
        electrode_potentials[:, i] = flatten_blocks(potential_block.(ele, centers,
                                                                     ranges, xorders))
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

function create_frame(sites)
    # TODO
    # Maybe add non-zero level offset between sites
    targets = [x2_block(site.x2, site.xorder) .+ yz_block(site.x2 / 2, site.xorder)
               for site in sites]

    zero_blocks = [zero_block(site.xorder) for site in sites]
    one_blocks = [c_block(1, site.xorder) for site in sites]
    single_site_line(idx, v) = [(i == idx ? v : zero_blocks[i])
                                for i in 1:length(sites)]

    freedom = Tuple{Vector{Vector{Float64}},Float64,Float64}[]
    dist = sites_gap(sites)
    for (i, site) in enumerate(sites)
        if site.level_min != 0 || site.level_max != 0
            push!(freedom, (single_site_line(i, one_blocks[i]), site.level_min,
                            site.level_max))
        end
        shape_relax = site.shape_relax
        for order in 6:2:site.xorder
            push!(freedom, (single_site_line(i, xn_block(order, 1, site.xorder)),
                            0, Inf))
        end
        if shape_relax == 0
            continue
        end
        x2 = site.x2
        yz = x2 / 2
        # Still assume single site is on the left
        if site.is_chain
            push!(freedom, (single_site_line(i, dx_block(1, site.xorder)),
                            -1250 * shape_relax, 1250 * shape_relax))
            push!(freedom, (single_site_line(i, zx_block(1, site.xorder)),
                            -yz * shape_relax / 4, yz * shape_relax / 4))
            push!(freedom, (single_site_line(i, z2_block(1, site.xorder)),
                            -yz * shape_relax / 4, yz * shape_relax / 4))
            push!(freedom, (single_site_line(i, x2_block(1, site.xorder)),
                            -x2 * shape_relax / 5, x2 * shape_relax / 5))
            push!(freedom, (single_site_line(i, xn_block(3, 1, site.xorder)),
                            -x2 * shape_relax * 0.2, x2 * shape_relax * 1))
            push!(freedom, (single_site_line(i, xn_block(4, 1, site.xorder)),
                            0, Inf))
        else
            push!(freedom, (single_site_line(i, dx_block(1, site.xorder)),
                            -1250 * shape_relax, 1250 * shape_relax))
            push!(freedom, (single_site_line(i, zx_block(1, site.xorder)),
                            -yz * shape_relax / 3, yz * shape_relax / 3))
            push!(freedom, (single_site_line(i, z2_block(1, site.xorder)),
                            -yz * shape_relax / 3, yz * shape_relax / 3))
            push!(freedom, (single_site_line(i, x2_block(1, site.xorder)),
                            -x2 * shape_relax / 8, x2 * shape_relax / 8))
            push!(freedom, (single_site_line(i, xn_block(3, 1, site.xorder)),
                            -x2 * shape_relax * 1, x2 * shape_relax * 0.3))
            push!(freedom, (single_site_line(i, xn_block(4, 1, site.xorder)),
                            dist > 350 ? -Inf : 0, Inf))
        end
    end
    push!(freedom, (one_blocks, -Inf, Inf))

    return construct_frame(get_electrodes(sites), targets,
                           freedom,
                           [(site.xidx_min, site.xidx_max) for site in sites],
                           [get_rf_center(site.xpos_um) for site in sites],
                           [site.xorder for site in sites])
end

const prefix = joinpath(@__DIR__, "../data/merge_coeff_20221021")

function dump_frame(sites, frame)
    neles = length(frame.electrodes)
    @assert size(frame.electrode_potentials, 2) == neles

    target = frame.target
    free_terms = Vector{Float64}[]
    free_term_limits = NTuple{2,Float64}[]
    for flex in frame.freedom
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
        push!(free_terms, frame.electrode_potentials \ flex[1])
        push!(free_term_limits, flex[2])
    end
    v = frame.electrode_potentials \ target
    coeff = frame.electrode_potentials
    B = qr(coeff').Q[:, (size(coeff, 1) + 1):end]
    nv, nt = size(B)
    @assert length(v) == nv

    return Dict("electrodes"=>frame.electrodes,
                "solution"=>v,
                "limited_solution"=>free_terms,
                "limits"=>[v[i] for v in free_term_limits, i in 1:2],
                "free_solution"=>B,
                "xpos_um"=>[site.xpos_um for site in sites])
end

function shuttling_frames()
    res = Dict{String,Any}[]
    for i in 0:shuttling_steps
        println("Shuttling: $i/$shuttling_steps")
        @time begin
            sites = sites_shuttling(i)
            populate_sites!(sites)
            frame = create_frame(sites)
            push!(res, dump_frame(sites, frame))
        end
    end
    return res
end

function merging_frames()
    res = Dict{String,Any}[]
    for i in 1:merging_steps
        println("Merging: $i/$merging_steps")
        @time begin
            sites = sites_merging(i)
            populate_sites!(sites)
            frame = create_frame(sites)
            push!(res, dump_frame(sites, frame))
        end
    end
    return res
end

function moveback_frames()
    res = Dict{String,Any}[]
    for i in 1:moveback_steps
        println("Moveback: $i/$moveback_steps")
        @time begin
            sites = sites_moveback(i)
            populate_sites!(sites)
            frame = create_frame(sites)
            push!(res, dump_frame(sites, frame))
        end
    end
    return res
end

const results = [shuttling_frames(), merging_frames(), moveback_frames()]

matopen("$(prefix).mat", "w") do mat
    write(mat, "data", results)
    write(mat, "electrode_names", solution.electrode_names)
end
