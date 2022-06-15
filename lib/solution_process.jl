#!/usr/bin/julia

module ProcessSolution

import ..PolyFit
using ..VoltageSolutions
import ..gradient
using ..OutputFiles
using NLsolve
using LinearAlgebra
using DelimitedFiles

struct FitCache{N,A<:AbstractArray{T,N} where T}
    fitter::PolyFit.PolyFitter{N}
    data::A
    cache::Array{PolyFit.PolyFitResult{N},N}
end

function FitCache(fitter::PolyFit.PolyFitter{N}, data::A) where (A<:AbstractArray{T,N} where T) where N
    cache = Array{PolyFit.PolyFitResult{N},N}(undef, size(data) .- fitter.orders)
    return FitCache{N,A}(fitter, data, cache)
end

function get_internal(cache::FitCache{N}, idx::NTuple{N,Integer}) where N
    if isassigned(cache.cache, idx...)
        return cache.cache[idx...]
    end
    data = @view cache.data[((:).(idx, idx .+ cache.fitter.orders))...]
    res = cache.fitter \ data
    cache.cache[idx...] = res
    return res
end

# pos is in index unit
function _best_fit_idx(ntotal, order, pos)
    # for fit at `index`, the data covered is `index:(index + order)`
    # with the center at index `index + order / 2`.
    # Therefore, the ideal index to use for `pos` is `pos - order / 2`
    idx = round(Int, pos - order / 2)
    if idx <= 1
        return 1
    elseif idx >= ntotal - order
        return ntotal - order
    end
    return idx
end

function get_shifted(cache::FitCache{N}, pos::NTuple{N}) where N
    sizes = size(cache.data)
    orders = cache.fitter.orders
    idxs = _best_fit_idx.(sizes, orders, pos)
    fit = get_internal(cache, idxs)
    return PolyFit.shift(fit, pos .- orders ./ 2 .- idxs)
end

function Base.get(cache::FitCache{N}, idx::NTuple{N}) where N
    return get_shifted(cache, idx)
end

function gradient(cache::FitCache{N}, dim, pos::Vararg{Any,N}) where N
    sizes = size(cache.data)
    orders = cache.fitter.orders
    idxs = _best_fit_idx.(sizes, orders, pos)
    fit = get_internal(cache, idxs)
    # center of the fit in index: idx + order / 2
    # position within fit: pos - idx - order / 2
    pos = pos .- idxs .- orders ./ 2
    return gradient(fit, dim, pos...)
end

function find_flat_point(data::A; init=ntuple(i->(size(data, i) + 1) / 2, Val(N))) where (A<:AbstractArray{T,N} where T) where N
    fitter = PolyFit.PolyFitter(ntuple(i->3, Val(N))...)
    cache = FitCache(fitter, data)
    function model!(g, x)
        xt = ntuple(i->x[i], Val(N))
        for i in 1:N
            g[i] = gradient(cache, i, xt...)
        end
    end
    res = nlsolve(model!, collect(init))
    return ntuple(i->res.zero[i], Val(N))
end

function find_all_flat_points(all_data::A; init=ntuple(i->(size(all_data, i) + 1) / 2, Val(N - 1))) where (A<:AbstractArray{T,N} where T) where N

    npoints = size(all_data, 3)
    all_res = Matrix{Float64}(undef, npoints, N - 1)

    for i in 1:npoints
        init = find_flat_point(@view(all_data[:, :, i]), init=init)
        all_res[i, :] .= init
    end
    return all_res
end

# X position of electrodes
# We need to know the axial (X) positions of the electrodes so that we can figure out
# which electrodes to use for generating potentials at a given location.

# From looking at the Sandia 3D model, each inner electrode is 70um wide in X direction
# (67um + 3um gap) and each outer quantum is 2x this (140um total).
# All the odd electrode are always located at the same X position as
# the electrode numbered one less than it.

# In unit of 70um and showing only the even electrodes, the order/positions
# of the electrodes are,

# Outer: |                47(O0)              |22(Q44-64)|          14(O0)          |
# Inner: |2(gap)|10(GND)|5(L0-8)|30(S10-0 x 5)|22(Q0-42) |8(S0-10,0-3)|6(GND)|2(gap)|

# where the number outside the parenthesis is the width in unit of 70um
# and the parenthesis marks the corresponding (even) electrode.
# The origin is located in the middle of the quantum region (11 from left and right).
# This distribution is cross-checked with the potential data
# by setting two pairs of electrode to opposite values and measuring the position
# of the zero crossing.

struct ElectrodePosition
    name::String
    left::Float64
    right::Float64
    up::Bool
end

function distance(pos::ElectrodePosition, x)
    if x < pos.left
        return pos.left - x
    elseif x > pos.right
        return x - pos.right
    else
        return 0.0
    end
end

const outer_positions = ElectrodePosition[]
const inner_positions = ElectrodePosition[]

function __populate_positions()
    begin_gnd = 12
    nL = 5
    nS = 6
    S_rep1 = 5
    nQ = 22
    S_rep2 = 1
    end_gnd = 8

    unit_um = 70

    @assert nQ % 2 == 0
    nQ_outer = nQ ÷ 2
    left_edge = -(begin_gnd + nL + nS * S_rep1 + nQ ÷ 2)

    pos_inner = left_edge + begin_gnd
    pos_outer = left_edge

    # Loading
    for i in 1:nL
        push!(inner_positions,
              ElectrodePosition("L$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(inner_positions,
              ElectrodePosition("L$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Transition 1
    for j in 1:S_rep1
        for i in nS:-1:1
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, true))
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, false))
            pos_inner += 1
        end
    end

    # Outer 1
    push!(outer_positions,
          ElectrodePosition("O0", pos_outer * unit_um,
                            pos_inner * unit_um, true))
    push!(outer_positions,
          ElectrodePosition("O1", pos_outer * unit_um,
                            pos_inner * unit_um, false))
    pos_outer = pos_inner

    # Quantum inner
    for i in 1:nQ
        push!(inner_positions,
              ElectrodePosition("Q$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(inner_positions,
              ElectrodePosition("Q$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Quantum outer
    for i in 1:nQ_outer
        i += nQ
        push!(outer_positions,
              ElectrodePosition("Q$(i * 2 - 2)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um, true))
        push!(outer_positions,
              ElectrodePosition("Q$(i * 2 - 1)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um, false))
        pos_outer += 2
    end
    @assert pos_inner == pos_outer

    # Transition 2
    for j in 1:S_rep2
        for i in 1:nS
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, true))
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, false))
            pos_inner += 1
        end
    end

    # S0-S3 appeared again at the end (shouldn't really matter......)
    for i in 1:4
        push!(inner_positions,
              ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(inner_positions,
              ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Outer 2
    push!(outer_positions,
          ElectrodePosition("O0", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um, true))
    push!(outer_positions,
          ElectrodePosition("O1", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um, false))
end

__populate_positions()

struct ConstraintSolution
    electrodes::Int
    nx::Int
    ny::Int
    nz::Int
    stride::NTuple{3,Float64}
    origin::NTuple{3,Float64}
    data::Array{Float64,4}
    electrode_index::Dict{String,Int}
    electrode_names::Vector{Vector{String}}
end

function ConstraintSolution(raw::VoltageSolution, aliases::Dict{Int,Int})
    origin_names = VoltageSolutions.electrode_names
    # Compute the mapping between id's
    id_map = zeros(Int, raw.electrodes)
    id = 0
    new_electrodes = raw.electrodes - length(aliases)
    data = Array{Float64}(undef, raw.nz, raw.ny, raw.nx, new_electrodes)
    electrode_names = Vector{Vector{String}}(undef, new_electrodes)
    electrode_index = Dict{String,Int}()
    for i in 1:raw.electrodes
        if i in keys(aliases)
            continue
        end
        id += 1
        id_map[i] = id
        data[:, :, :, id] .= @view raw.data[:, :, :, i]
        name = origin_names[i]
        electrode_index[name] = id
        electrode_names[id] = [name]
    end
    @assert new_electrodes == id
    for (k, v) in aliases
        # The user should connect directly to the final one
        @assert !(v in keys(aliases))
        id = id_map[v]
        @assert id != 0
        data[:, :, :, id] .= @view(data[:, :, :, id]) .+ @view(raw.data[:, :, :, k])
        name = origin_names[k]
        electrode_index[name] = id
        push!(electrode_names[id], name)
    end
    return ConstraintSolution(new_electrodes, raw.nx, raw.ny, raw.nz,
                              raw.stride, raw.origin,
                              data, electrode_index, electrode_names)
end

function ConstraintSolution(raw::VoltageSolution, aliases::Dict{String,String})
    mapping = VoltageSolutions.electrode_index
    return ConstraintSolution(raw, Dict(mapping[k]=>mapping[v] for (k, v) in aliases))
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        export $(Symbol(name, "_index_to_axis"))
        VoltageSolutions.$(Symbol(name, "_index_to_axis"))(sol::ConstraintSolution, i) = (i - 1) * sol.stride[$i] + sol.origin[$i]
        export $(Symbol(name, "_axis_to_index"))
        VoltageSolutions.$(Symbol(name, "_axis_to_index"))(sol::ConstraintSolution, a) = (a - sol.origin[$i]) / sol.stride[$i] + 1
    end
end

mutable struct ElectrodeSearchState
    inner_idx1::Int
    inner_idx2::Int
    outer_idx1::Int
    outer_idx2::Int

    const pos::Float64
    const inner_candidates::Vector{ElectrodePosition}
    const outer_candidates::Vector{ElectrodePosition}
    function ElectrodeSearchState(pos)
        inner_idx2 = searchsortedfirst(inner_positions, pos, lt=(x, y)->x.right < y)
        outer_idx2 = searchsortedfirst(outer_positions, pos, lt=(x, y)->x.right < y)

        return new(inner_idx2 - 1, inner_idx2, outer_idx2 - 1, outer_idx2,
                   pos, ElectrodePosition[], ElectrodePosition[])
    end
end

function _find_next_distance!(state::ElectrodeSearchState)
    # First find the closest distance
    dist = Inf
    if state.inner_idx2 <= length(inner_positions)
        dist = min(dist, distance(inner_positions[state.inner_idx2], state.pos))
    end
    if state.inner_idx1 > 0
        dist = min(dist, distance(inner_positions[state.inner_idx1], state.pos))
    end
    if state.outer_idx2 <= length(outer_positions)
        dist = min(dist, distance(outer_positions[state.outer_idx2], state.pos))
    end
    if state.outer_idx1 > 0
        dist = min(dist, distance(outer_positions[state.outer_idx1], state.pos))
    end
    if !isfinite(dist)
        error("Unable to find enough terms")
    end
    empty!(state.inner_candidates)
    empty!(state.outer_candidates)
    while state.inner_idx2 <= length(inner_positions)
        epos = inner_positions[state.inner_idx2]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.inner_candidates, epos)
        state.inner_idx2 += 1
    end
    while state.inner_idx1 > 0
        epos = inner_positions[state.inner_idx1]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.inner_candidates, epos)
        state.inner_idx1 -= 1
    end
    while state.outer_idx2 <= length(outer_positions)
        epos = outer_positions[state.outer_idx2]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.outer_candidates, epos)
        state.outer_idx2 += 1
    end
    while state.outer_idx1 > 0
        epos = outer_positions[state.outer_idx1]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.outer_candidates, epos)
        state.outer_idx1 -= 1
    end
    @assert !isempty(state.inner_candidates) || !isempty(state.outer_candidates)
end

"""
Find exactly n electrodes that are the closest in axial (X) position to `pos`.

The constraints from `solution` will be used to make sure
"""
function find_n_electrodes(solution::ConstraintSolution, pos, n)
    nup = 0
    ndown = 0
    res = Set{Int}()

    search_state = ElectrodeSearchState(pos)

    inner_up = Int[]
    inner_down = Int[]
    outer_up = Int[]
    outer_down = Int[]
    while true
        nleft = n - length(res)
        if nleft <= 0
            return res
        end

        _find_next_distance!(search_state)

        empty!(inner_up)
        empty!(inner_down)
        for p in search_state.inner_candidates
            id = solution.electrode_index[p.name]
            # Ground or already added
            if id == 1
                continue
            elseif id in res
                # Don't add duplicate electrode but do count up vs down
                if p.up
                    nup += 1
                else
                    ndown += 1
                end
                continue
            end
            push!(p.up ? inner_up : inner_down, id)
        end
        n_inner_up = length(inner_up)
        n_inner_down = length(inner_down)
        if n_inner_up + n_inner_down <= nleft
            nleft -= n_inner_up + n_inner_down
            nup += n_inner_up
            ndown += n_inner_down
            union!(res, inner_up)
            union!(res, inner_down)
        else
            # Found too many: try to balance up and down
            # Do note that nup + ndown could be more than n due to shorted electrodes
            ndown_use = min(max((nleft + nup - ndown) ÷ 2, 0), n_inner_down, nleft)
            nup_use = nleft - ndown_use
            union!(res, @view inner_up[1:nup_use])
            union!(res, @view inner_down[1:ndown_use])
            return res
        end

        empty!(outer_up)
        empty!(outer_down)
        for p in search_state.outer_candidates
            id = solution.electrode_index[p.name]
            # Ground or already added
            if id == 1
                continue
            elseif id in res
                # Don't add duplicate electrode but do count up vs down
                if p.up
                    nup += 1
                else
                    ndown += 1
                end
                continue
            end
            push!(p.up ? outer_up : outer_down, id)
        end
        n_outer_up = length(outer_up)
        n_outer_down = length(outer_down)
        if n_outer_up + n_outer_down <= nleft
            nleft -= n_outer_up + n_outer_down
            nup += n_outer_up
            ndown += n_outer_down
            union!(res, outer_up)
            union!(res, outer_down)
        else
            # Found too many: try to balance up and down
            # Do note that nup + ndown could be more than n due to shorted electrodes
            ndown_use = min(max((nleft + nup - ndown) ÷ 2, 0), n_outer_down, nleft)
            nup_use = nleft - ndown_use
            union!(res, @view outer_up[1:nup_use])
            union!(res, @view outer_down[1:ndown_use])
            return res
        end
    end
end

const _subarray_T = typeof(@view zeros(0, 0, 0, 1)[:, :, :, 1])

struct ElectrodesFitCache
    fitter::PolyFit.PolyFitter{3}
    solution::ConstraintSolution
    cache::Vector{FitCache{3,_subarray_T}}
    function ElectrodesFitCache(fitter::PolyFit.PolyFitter{3},
                                solution::ConstraintSolution)
        return new(fitter, solution,
                   Vector{FitCache{3,_subarray_T}}(undef, solution.electrodes))
    end
end

function Base.get(cache::ElectrodesFitCache, idx)
    if isassigned(cache.cache, idx)
        return cache.cache[idx]
    end
    fit_cache = FitCache(cache.fitter, @view cache.solution.data[:, :, :, idx])
    cache.cache[idx] = fit_cache
    return fit_cache
end

Base.get(cache::ElectrodesFitCache, name::AbstractString) =
    get(cache, cache.solution.electrode_index[name])

Base.get(cache::ElectrodesFitCache, electrode, pos::NTuple{3}) =
    get(get(cache, electrode), pos)

# Terms we care about
# x, y, z, xy, yz, xz, z^2 - y^2, x^2 - (y^2 + z^2) / 2, x^3, x^4
# Since we care about the symmetry of the x^2 and z^2 term,
# we actually do need to scale the x, y and z correctly.
# stride should be in um, voltage should be in V
function get_compensate_terms1(res::PolyFit.PolyFitResult{3}, stride)
    # axis order of fitting result is z, y, x
    # axis order of stride is x, y, z
    raw_x = res[0, 0, 1]
    raw_y = res[0, 1, 0]
    raw_z = res[1, 0, 0]

    raw_xy = res[0, 1, 1]
    raw_yz = res[1, 1, 0]
    raw_zx = res[1, 0, 1]

    raw_x2 = res[0, 0, 2]
    raw_y2 = res[0, 2, 0]
    raw_z2 = res[2, 0, 0]

    raw_x3 = res[0, 0, 3]
    raw_x4 = res[0, 0, 4]

    scaled_x = raw_x / stride[1]
    scaled_y = raw_y / stride[2]
    scaled_z = raw_z / stride[3]

    scaled_xy = raw_xy / stride[1] / stride[2]
    scaled_yz = raw_yz / stride[2] / stride[3]
    scaled_zx = raw_zx / stride[3] / stride[1]

    scaled_x2 = raw_x2 / stride[1]^2
    scaled_y2 = raw_y2 / stride[2]^2
    scaled_z2 = raw_z2 / stride[3]^2

    scaled_x3 = raw_x3 / stride[1]^3
    scaled_x4 = raw_x4 / stride[1]^4

    xx = (scaled_x2 - scaled_y2 - scaled_z2) / 3
    zz = (scaled_z2 - scaled_y2) / 2

    # Current units are V/um^n
    # Expected units
    # DX/DY/DZ: V/m
    # XY, YZ, ZX, ZZ, XX: 525 uV / (2.74 um)^2
    # X3: 525 uV / (2.74 um)^3
    # X4: 525 uV / (2.74 um)^4
    scale_1 = 1e6
    scale_2 = (2.74^2 / 525e-6)
    scale_3 = (2.74^3 / 525e-6)
    scale_4 = (2.74^4 / 525e-6)
    return (dx=scaled_x * scale_1, dy=scaled_y * scale_1, dz=scaled_z * scale_1,
            xy=scaled_xy * scale_2, yz=scaled_yz * scale_2, zx=scaled_zx * scale_2,
            z2=zz * scale_2, x2=xx * scale_2, x3=scaled_x3 * scale_3,
            x4=scaled_x4 * scale_4)
end

function solve_terms1(fits::Vector{PolyFit.PolyFitResult{3}}, stride)
    nfits = length(fits)
    coefficient = Matrix{Float64}(undef, 10, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(get_compensate_terms1(fits[i], stride))
    end
    X = coefficient \ Matrix(I, 10, 10)
    @assert size(X, 2) == 10
    return (dx=X[:, 1], dy=X[:, 2], dz=X[:, 3],
            xy=X[:, 4], yz=X[:, 5], zx=X[:, 6],
            z2=X[:, 7], x2=X[:, 8], x3=X[:, 9], x4=X[:, 10])
end

function compensate_fitter1(solution::ConstraintSolution)
    fitter = PolyFit.PolyFitter(2, 2, 4)
    return ElectrodesFitCache(fitter, solution)
end

function get_compensate_terms1(cache::ElectrodesFitCache, pos::NTuple{3})
    # pos is in xyz index

    x_coord = x_index_to_axis(cache.solution, pos[1]) .* 1000
    ele_select = find_n_electrodes(cache.solution, x_coord, 20)
    ele_select = sort!(collect(ele_select))
    fits = [get(cache, e, (pos[3], pos[2], pos[1])) for e in ele_select]
    # Change stride to um in unit
    return ele_select, solve_terms1(fits, cache.solution.stride .* 1000)
end

struct CenterTracker
    zy_index::Matrix{Float64}
end

function Base.get(tracker::CenterTracker, xidx)
    # return (y, z)
    nx = size(tracker.zy_index, 1)
    lb_idx = min(max(floor(Int, xidx), 1), nx)
    ub_idx = min(max(ceil(Int, xidx), 1), nx)
    y_lb = tracker.zy_index[lb_idx, 2]
    z_lb = tracker.zy_index[lb_idx, 1]
    if lb_idx == ub_idx
        return y_lb, z_lb
    end
    @assert ub_idx == lb_idx + 1
    y_ub = tracker.zy_index[ub_idx, 2]
    z_ub = tracker.zy_index[ub_idx, 1]
    c_ub = xidx - lb_idx
    c_lb = ub_idx - xidx
    return y_lb * c_lb + y_ub * c_ub, z_lb * c_lb + z_ub * c_ub
end

function load_short_map(fname)
    m = readdlm(fname, ',', String)
    res = Dict{String,String}()
    for i in 1:size(m, 1)
        res[m[i, 1]] = m[i, 2]
    end
    return res
end

function _get_data_line(solution::ConstraintSolution, mapfile::MapFile,
                        electrodes, term)
    term_map = Dict(zip(electrodes, term))
    nelectrodes = length(mapfile.names)
    values = Vector{Float64}(undef, nelectrodes)
    for i in 1:nelectrodes
        id = get(solution.electrode_index, mapfile.names[i], -1)
        values[i] = get(term_map, id, 0.0)
    end
    return values
end

function compensation_to_file(solution::ConstraintSolution, mapfile::MapFile,
                              electrodes, terms)
    term_names = String[]
    term_values = Vector{Float64}[]
    nelectrodes = length(mapfile.names)
    # DX
    push!(term_names, "DX")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.dx))
    # DY
    push!(term_names, "DY")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.dy))
    # DZ
    push!(term_names, "DZ")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.dz))
    # QZY
    push!(term_names, "QZY")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.yz))
    # QZZ
    push!(term_names, "QZZ")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.z2))
    # QXZ
    push!(term_names, "QXZ")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.zx))
    # X1
    # DX is in V/m, X1 is in 525 uV / 2.74 um
    push!(term_names, "X1")
    push!(term_values, _get_data_line(solution, mapfile, electrodes,
                                      terms.dx .* (525 / 2.74)))
    # X2
    push!(term_names, "X2")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.x2))
    # X3
    push!(term_names, "X3")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.x3))
    # X4
    push!(term_names, "X4")
    push!(term_values, _get_data_line(solution, mapfile, electrodes, terms.x4))
    # JunctionCenter
    push!(term_names, "JunctionCenter")
    push!(term_values, zeros(nelectrodes))
    # JunctionTransition
    push!(term_names, "JunctionTransition")
    push!(term_values, zeros(nelectrodes))
    # LoadInners
    push!(term_names, "LoadInners")
    push!(term_values, zeros(nelectrodes))
    # LoadOuters
    push!(term_names, "LoadOuters")
    push!(term_values, zeros(nelectrodes))
    # Eject
    push!(term_names, "Eject")
    push!(term_values, zeros(nelectrodes))

    return CompensationFile(mapfile, term_names, term_values)
end

end
