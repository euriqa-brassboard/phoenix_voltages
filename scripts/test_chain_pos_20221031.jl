#!/usr/bin/julia

using JuMP
using Interpolations

mutable struct ModelBuilder{F,∇F,∇²F}
    model::Model
    f::F # Potential function
    ∇f::∇F
    ∇²f::∇²F
    ions::Vector{Pair{VariableRef,Float64}}
    last_ion::VariableRef
    last_barrier::Float64
    barrier_after_ion::Bool
    function ModelBuilder{F,∇F,∇²F}(model::Model, f::F, ∇f::∇F, ∇²f::∇²F) where {F,∇F,∇²F}
        builder = new{F,∇F,∇²F}(model, f, ∇f, ∇²f, Pair{VariableRef,Float64}[])
        register(model, :potential, 1, f, ∇f, ∇²f)
        builder.barrier_after_ion = false
        return builder
    end
    function ModelBuilder(model::Model, f::F, ∇f::∇F, ∇²f::∇²F) where {F,∇F,∇²F}
        return ModelBuilder{F,∇F,∇²F}(model, f, ∇f, ∇²f)
    end
end

function add_barrier!(builder::ModelBuilder, barrier)
    model = builder.model
    builder.last_barrier = barrier
    if !builder.barrier_after_ion
        if isdefined(builder, :last_ion)
            set_upper_bound(builder.last_ion, barrier)
        end
        builder.barrier_after_ion = true
    end
    return
end

function add_ion!(builder::ModelBuilder, init_value, charge)
    model = builder.model
    ion = @variable(model)
    set_start_value(ion, init_value)
    push!(builder.ions, ion=>charge)
    if builder.barrier_after_ion
        set_lower_bound(ion, builder.last_barrier)
        builder.barrier_after_ion = false
    elseif isdefined(builder, :last_ion)
        @constraint(model, builder.last_ion <= ion)
    end
    builder.last_ion = ion
    return ion
end

function finalize_model!(builder::ModelBuilder)
    model = builder.model
    obj = @NLexpression(model, 0)
    nions = length(builder.ions)
    for (i1, (ion1, charge1)) in enumerate(builder.ions)
        obj = @NLexpression(model, obj + potential(ion1))
        for i2 in (i1 + 1):nions
            ion2, charge2 = builder.ions[i2]
            obj = @NLexpression(model, obj + charge1 * charge2 / (ion2 - ion1))
        end
    end
    @NLobjective(model, Min, obj)
    return model
end

function interpolate_ref_functions(fi::Ref)
    return (x->fi[](x), x->only(Interpolations.gradient(fi[], x)),
            x->only(Interpolations.hessian(fi[], x)))
end
