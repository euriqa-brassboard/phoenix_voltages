module Optimizers

import NLopt
using JuMP
using LinearAlgebra

function gen_minmax_model(B, x0)
    nx, nt = size(B)
    @assert nx == length(x0)

    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    @variable(model, t[1:nt])
    @variable(model, maxv)
    @constraint(model, maxv .>= B * t .+ x0)
    @constraint(model, maxv .>= .-(B * t .+ x0))
    @objective(model, Min, maxv)
    return model, t
end

function optimize_minmax_span(B, x0)
    model, t = gen_minmax_model(B, x0)
    JuMP.optimize!(model)
    return B * value.(t) .+ x0
end

function optimize_minmax(A, y::AbstractVector)
    x0 = A \ y
    ny, nx = size(A)
    nt = nx - ny
    if nt <= 0
        return x0
    end
    # With the A * x = y constraints,
    # the degrees of freedom left in x are the ones that satisfies A * x = 0
    # In another word, these are the x's that are orthogonal to all rows of A.
    # We can find the basis set that spans such space using QR decomposition.
    B = qr(A').Q[:, (ny + 1):nx]
    model, t = gen_minmax_model(B, x0)
    JuMP.optimize!(model)
    return B * value.(t) .+ x0
end

function optimize_minmax(A, y::AbstractMatrix)
    ny, nx = size(A)
    @assert ny == size(y, 1)
    ns = size(y, 2)

    x = Matrix{Float64}(undef, nx, ns)
    for i in 1:ns
        x[:, i] .= optimize_minmax(A, @view(y[:, i]))
    end
    return x
end

end
