#!/usr/bin/julia

using JuMP
using Interpolations

struct IonChainModel{F,∇F,∇²F}
    model::Model
    f::F # Potential function
    ∇f::∇F
    ∇²f::∇²F
    ions::Vector{Pair{VariableRef,Float64}}
    function IonChainModel{F,∇F,∇²F}(model::Model, f::F, ∇f::∇F, ∇²f::∇²F) where {F,∇F,∇²F}
        register(model, :potential, 1, f, ∇f, ∇²f)
        return new{F,∇F,∇²F}(model, f, ∇f, ∇²f, Pair{VariableRef,Float64}[])
    end
    function IonChainModel(model::Model, f::F, ∇f::∇F, ∇²f::∇²F) where {F,∇F,∇²F}
        return IonChainModel{F,∇F,∇²F}(model, f, ∇f, ∇²f)
    end
end

mutable struct ModelBuilder{F,∇F,∇²F}
    model::IonChainModel{F,∇F,∇²F}
    barrier_after_ion::Bool
    last_ion::VariableRef
    last_barrier::Float64
    function ModelBuilder{F,∇F,∇²F}(model::IonChainModel{F,∇F,∇²F}) where {F,∇F,∇²F}
        return new{F,∇F,∇²F}(model, false)
    end
    function ModelBuilder(model::IonChainModel{F,∇F,∇²F}) where {F,∇F,∇²F}
        return ModelBuilder{F,∇F,∇²F}(model)
    end
end

function add_barrier!(builder::ModelBuilder, barrier)
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
    m = model.model
    ion = @variable(m)
    set_start_value(ion, init_value)
    push!(model.ions, ion=>charge)
    if builder.barrier_after_ion
        set_lower_bound(ion, builder.last_barrier)
        builder.barrier_after_ion = false
    elseif isdefined(builder, :last_ion)
        @constraint(m, builder.last_ion <= ion)
    end
    builder.last_ion = ion
    return ion
end

function finalize_model!(builder::ModelBuilder)
    model = builder.model
    m = model.model
    obj = @NLexpression(m, 0)
    ions = model.ions
    nions = length(ions)
    for (i1, (ion1, charge1)) in enumerate(ions)
        obj = @NLexpression(m, obj + potential(ion1))
        for i2 in (i1 + 1):nions
            ion2, charge2 = ions[i2]
            obj = @NLexpression(m, obj + charge1 * charge2 / (ion2 - ion1))
        end
    end
    @NLobjective(m, Min, obj)
    return
end

function update_init_pos!(model::IonChainModel)
    # Setting start value clears the value
    values = [value(ion) for (ion, charge) in model.ions]
    for ((ion, charge), v) in zip(model.ions, values)
        set_start_value(ion, v)
    end
end

function interpolate_ref_functions(fi::Ref)
    return (x->fi[](x), x->only(Interpolations.gradient(fi[], x)),
            x->only(Interpolations.hessian(fi[], x)))
end
