#!/usr/bin/julia

module ProcessSolution

import ..PolyFit
import ..gradient
using NLsolve

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

function Base.get(cache::FitCache{N}, idx::NTuple{N,Integer}) where N
    if checkbounds(Bool, cache.cache, idx...)
        return get_internal(cache, idx)
    end
    return get_shifted(cache, idx)
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

# Outer: |                47(O0)              |22(Q44-64)|        14(O0)        |
# Inner: |2(gap)|10(GND)|5(L0-8)|30(S10-0 x 5)|22(Q0-42) |6(S0-10)|6(GND)|2(gap)|

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
                                (pos_inner + 1) * unit_um))
        push!(inner_positions,
              ElectrodePosition("L$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um))
        pos_inner += 1
    end

    # Transition 1
    for j in 1:S_rep1
        for i in nS:-1:1
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um))
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um))
            pos_inner += 1
        end
    end

    # Outer 1
    push!(outer_positions,
          ElectrodePosition("O0", pos_outer * unit_um,
                            pos_inner * unit_um))
    push!(outer_positions,
          ElectrodePosition("O1", pos_outer * unit_um,
                            pos_inner * unit_um))
    pos_outer = pos_inner

    # Quantum inner
    for i in 1:nQ
        push!(inner_positions,
              ElectrodePosition("Q$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um))
        push!(inner_positions,
              ElectrodePosition("Q$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um))
        pos_inner += 1
    end

    # Quantum outer
    for i in 1:nQ_outer
        i += nQ
        push!(outer_positions,
              ElectrodePosition("Q$(i * 2 - 2)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um))
        push!(outer_positions,
              ElectrodePosition("Q$(i * 2 - 1)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um))
        pos_outer += 2
    end
    @assert pos_inner == pos_outer

    # Transition 2
    for j in 1:S_rep2
        for i in 1:nS
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um))
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um))
            pos_inner += 1
        end
    end

    # Outer 1
    push!(outer_positions,
          ElectrodePosition("O0", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um))
    push!(outer_positions,
          ElectrodePosition("O1", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um))
end

__populate_positions()

end
