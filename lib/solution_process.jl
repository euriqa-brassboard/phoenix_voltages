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
    return PolyFit.shift(fit, pos .- idxs)
end

function Base.get(cache::FitCache{N}, idx::NTuple{N,Integer}) where N
    if checkbounds(Bool, cache.cache, idx...)
        return get_internal(cache, idx)
    end
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

end
