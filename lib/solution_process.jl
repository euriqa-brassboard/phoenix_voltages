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

function Base.get(cache::FitCache{N}, idx::NTuple{N}) where N
    if isassigned(cache.cache, idx...)
        return cache.cache[idx...]
    end
    data = @view cache.data[((:).(idx, idx .+ cache.fitter.orders))...]
    res = cache.fitter \ data
    cache.cache[idx...] = res
    return res
end

# pos = 0 in the middle of the data range
function _best_fit_idx(ntotal, order, pos)
    # for fit at `index`, the data covered is `index:(index + order)`
    # with the center at index `index + order / 2`.
    # this corresponds to a position of `index + order / 2 - (ntotal + 1) / 2`
    # or `index - (ntotal + 1 - order) / 2`
    # Therefore, the ideal index to use for `pos` is `pos + (ntotal + 1 - order) / 2`
    idx = round(Int, pos + (ntotal + 1 - order) / 2)
    if idx <= 1
        return 1
    elseif idx >= ntotal - order
        return ntotal - order
    end
    return idx
end

function gradient(cache::FitCache{N}, dim, pos::Vararg{Any,N}) where N
    sizes = size(cache.data)
    orders = cache.fitter.orders
    idxs = _best_fit_idx.(sizes, orders, pos)
    fit = get(cache, idxs)
    # center of the fit in index: idx + order / 2
    # center of the fit in position: idx - (size + 1 - order) / 2
    # position within fit: pos - idx + (size + 1 - order) / 2
    pos = pos .- idxs .+ (sizes .+ 1 .- orders) ./ 2
    return gradient(fit, dim, pos...)
end

function find_flat_point(data::A; init=ntuple(i->0.0, Val(N))) where (A<:AbstractArray{T,N} where T) where N
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

end
