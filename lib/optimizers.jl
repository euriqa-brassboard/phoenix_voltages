module Optimizers

import NLopt
using JuMP
using LinearAlgebra

struct LinearConstraintMinMaxModel
    model::Model
    B::Matrix{VariableRef}
    t::Vector{VariableRef}
    x0::Vector{VariableRef}
end

function minmax_objective(x...)
    return max(abs.(x)...)
end

function gen_minmax_model_template(nx, nt)
    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    @variable(model, x0[1:nx])
    @variable(model, t[1:nt])
    @variable(model, B[1:nx, 1:nt])
    @variable(model, maxv)
    @constraint(model, maxv .>= B * t .+ x0)
    @constraint(model, maxv .>= .-(B * t .+ x0))
    @objective(model, Min, maxv)
    return LinearConstraintMinMaxModel(model, B, t, x0)
end

const model_cache = Dict{Tuple{Int,Int},LinearConstraintMinMaxModel}()

function gen_fixed_minmax_model(B, x0)
    nx, nt = size(B)
    @assert nx == length(x0)

    model = pop!(model_cache, (nx, nt), nothing)
    if model === nothing
        model = gen_minmax_model_template(nx, nt)
    end
    for i in 1:length(B)
        fix(model.B[i], B[i])
    end
    for i in 1:length(x0)
        fix(model.x0[i], x0[i])
    end
    return model
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
    model = gen_fixed_minmax_model(B, x0)
    JuMP.optimize!(model.model)
    t = Float64[value(t) for t in model.t]
    ny, nx = size(A)
    model_cache[(nx, ny)] = model
    return B * t .+ x0
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
