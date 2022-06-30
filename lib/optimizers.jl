module Optimizers

using NLopt
using JuMP

struct LinearConstraintMinMaxModel
    model::Model
    A::Matrix{VariableRef}
    x::Vector{VariableRef}
    y::Vector{VariableRef}
end

function minmax_objective(x...)
    return max(abs.(x)...)
end

function gen_minmax_model_template(nx, ny)
    model = Model(NLopt.Optimizer)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    @variable(model, x[1:nx])
    @variable(model, y[1:ny])
    @variable(model, A[1:ny, 1:nx])
    @constraint(model, A * x .== y)
    register(model, :minmax_objective, nx, minmax_objective, autodiff=true)
    @NLobjective(model, Min, minmax_objective(x...))
    return LinearConstraintMinMaxModel(model, A, x, y)
end

const model_cache = Dict{Tuple{Int,Int},LinearConstraintMinMaxModel}()

function gen_fixed_minmax_model(A, y)
    ny, nx = size(A)
    @assert ny == length(y)
    x = A \ y

    model = pop!(model_cache, (nx, ny), nothing)
    if model === nothing
        model = gen_minmax_model_template(nx, ny)
    end
    for i in 1:length(A)
        fix(model.A[i], A[i])
    end
    for i in 1:length(y)
        fix(model.y[i], y[i])
    end
    for i in 1:length(x)
        set_start_value(model.x[i], x[i])
    end
    return model
end

function optimize_minmax(A, y::AbstractVector)
    model = gen_fixed_minmax_model(A, y)
    JuMP.optimize!(model.model)
    res = Float64[value(x) for x in model.x]
    ny, nx = size(A)
    model_cache[(nx, ny)] = model
    return res
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
