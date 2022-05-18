#!/usr/bin/julia

# N dimensional multipole fitting
struct PolyFitter{N}
    orders::NTuple{N,Int}
    coefficient::Matrix{Float64}
    function PolyFitter(orders::Vararg{Integer,N}) where N
        sizes = orders .+ 1
        nterms = prod(sizes)
        coefficient = Matrix{Float64}(undef, nterms, nterms)
        lindices = LinearIndices(sizes)
        cindices = CartesianIndices(sizes)
        # Index for position
        for ipos in lindices
            # Position of the point, with the origin in the middle of the grid.
            pos = Tuple(cindices[ipos]) .- 1 .- orders ./ 2
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

function Base.:\(fitter::PolyFitter{N}, data::AbstractMatrix) where N
    return PolyFitResult{N}(fitter.orders, fitter.coefficient \ vec(data))
end

function Base.getindex(res::PolyFitResult{N}, order::Vararg{Integer,N}) where N
    i = LinearIndices(res.orders .+ 1)[CartesianIndex(order .+ 1)]
    return res.coefficient[i]
end

function Base.setindex!(res::PolyFitResult{N}, val, order::Vararg{Integer,N}) where N
    i = LinearIndices(res.orders .+ 1)[CartesianIndex(order .+ 1)]
    res.coefficient[i] = val
    return
end
