#!/usr/bin/julia

# N dimensional multipole fitting
struct PolyFitter{N}
    orders::NTuple{N,Int}
    coefficient::Matrix{Float64}
    function PolyFitter(orders::Vararg{Integer,N}; center=orders ./ 2) where N
        sizes = orders .+ 1
        nterms = prod(sizes)
        coefficient = Matrix{Float64}(undef, nterms, nterms)
        lindices = LinearIndices(sizes)
        cindices = CartesianIndices(sizes)
        # Index for position
        for ipos in lindices
            # Position of the point, with the origin in the middle of the grid.
            pos = Tuple(cindices[ipos]) .- 1 .- center
            # Index for the polynomial order
            for iorder in lindices
                order = Tuple(cindices[iorder]) .- 1
                coefficient[ipos, iorder] = prod(pos.^order)
            end
        end
        return new{N}(orders, coefficient)
    end
end

struct PolyFitResult{N}
    orders::NTuple{N,Int}
    coefficient::Vector{Float64}
end

function (res::PolyFitResult{N})(pos::Vararg{Any,N}) where N
    sizes = res.orders .+ 1
    lindices = LinearIndices(sizes)
    cindices = CartesianIndices(sizes)
    v = 0.0
    for iorder in lindices
        order = Tuple(cindices[iorder]) .- 1
        v += res.coefficient[iorder] * prod(pos.^order)
    end
    return v
end

# Gradient along dimension
function grad(res::PolyFitResult{N}, dim, pos::Vararg{Any,N}) where N
    sizes = res.orders .+ 1
    lindices = LinearIndices(sizes)
    cindices = CartesianIndices(sizes)
    order_offset = ntuple((i)->i == dim ? -1 : 0, Val(N))
    v = 0.0
    for iorder in lindices
        # Original polymomial orders along each dimensions
        order = Tuple(cindices[iorder]) .- 1
        factor = order[dim]
        if factor == 0
            continue
        end
        term = factor * res.coefficient[iorder]
        @inbounds for i in 1:N
            o = i == dim ? order[i] - 1 : order[i]
            term *= pos[i]^o
        end
        v += term
    end
    return v
end

function order_index(orders::NTuple{N,Integer}, order::Vararg{Integer,N}) where N
    return LinearIndices(orders .+ 1)[CartesianIndex(order .+ 1)]
end

function order_index(res::PolyFitResult{N}, order::Vararg{Integer,N}) where N
    return order_index(res.orders, order...)
end

function Base.:\(fitter::PolyFitter{N}, data::AbstractMatrix) where N
    return PolyFitResult{N}(fitter.orders, fitter.coefficient \ vec(data))
end

function Base.getindex(res::PolyFitResult{N}, order::Vararg{Integer,N}) where N
    return res.coefficient[order_index(res, order...)]
end

function Base.setindex!(res::PolyFitResult{N}, val, order::Vararg{Integer,N}) where N
    res.coefficient[order_index(res, order...)] = val
    return
end

shifted_term(max_order, term_order, shift) =
    shift^(max_order - term_order) * binomial(max_order, term_order)

function shifted_coefficient(res::PolyFitResult{N}, shift::NTuple{N},
                             order::Vararg{Integer,N}) where N
    v = 0.0
    sizes = res.orders .+ 1
    lindices = LinearIndices(sizes)
    cindices = CartesianIndices(sizes)
    for lidx in lindices
        cidx = cindices[lidx]
        term_order = Tuple(cidx) .- 1
        if !all(term_order .>= order)
            continue
        end
        v += prod(shifted_term.(term_order, order, shift)) * res.coefficient[lidx]
    end
    return v
end

function shift(res::PolyFitResult{N}, shift::NTuple{N}) where N
    coefficient = similar(res.coefficient)
    sizes = res.orders .+ 1
    lindices = LinearIndices(sizes)
    cindices = CartesianIndices(sizes)
    for lidx in lindices
        cidx = cindices[lidx]
        order = Tuple(cidx) .- 1
        coefficient[lidx] = shifted_coefficient(res, shift, order...)
    end
    return PolyFitResult{N}(res.orders, coefficient)
end
