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
