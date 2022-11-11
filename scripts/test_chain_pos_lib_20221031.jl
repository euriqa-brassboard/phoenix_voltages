#!/usr/bin/julia

using JuMP
using Interpolations
using LinearAlgebra

# TODO: support different ion experiencing different potentials
# (ions with different mass experiences different pseudo potential)

struct IonInfo
    pos::VariableRef
    charge::Float64
    mass::Float64
    IonInfo(pos, charge=1, mass=1) = new(pos, charge, mass)
end

struct IonChainModel{F,∇F,∇²F}
    model::Model
    f::F # Potential function
    ∇f::∇F
    ∇²f::∇²F
    ions::Vector{IonInfo}
    function IonChainModel{F,∇F,∇²F}(model::Model, f::F, ∇f::∇F, ∇²f::∇²F) where {F,∇F,∇²F}
        register(model, :potential, 1, f, ∇f, ∇²f)
        return new{F,∇F,∇²F}(model, f, ∇f, ∇²f, IonInfo[])
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

function add_ion!(builder::ModelBuilder, init_value, charge=1, mass=1)
    model = builder.model
    m = model.model
    ion = @variable(m)
    set_start_value(ion, init_value)
    push!(model.ions, IonInfo(ion, charge, mass))
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
    for (i1, ion1) in enumerate(ions)
        obj = @NLexpression(m, obj + potential(ion1.pos))
        for i2 in (i1 + 1):nions
            ion2 = ions[i2]
            obj = @NLexpression(m, obj + ion1.charge * ion2.charge / (ion2.pos - ion1.pos))
        end
    end
    @NLobjective(m, Min, obj)
    return
end

function update_init_pos!(model::IonChainModel)
    # Setting start value clears the value
    values = [value(ion.pos) for ion in model.ions]
    for (ion, v) in zip(model.ions, values)
        set_start_value(ion.pos, v)
    end
end

function axial_modes(model::IonChainModel)
    values = [value(ion.pos) for ion in model.ions]
    nions = length(values)
    H = zeros(nions, nions)
    for i in 1:nions
        pos = values[i]
        H[i, i] = model.∇²f(pos) / model.ions[i].mass
    end
    for i2 in 2:nions
        pos2 = values[i2]
        mass2 = model.ions[i2].mass
        for i1 in 1:i2 - 1
            pos1 = values[i1]
            mass1 = model.ions[i1].mass
            mass12 = sqrt(mass1 * mass2)
            term = 2 / (pos2 - pos1)^3
            H[i1, i1] += term / mass1
            H[i2, i2] += term / mass2
            H[i1, i2] -= term / mass12
            H[i2, i1] -= term / mass12
        end
    end
    evs = eigvals(H)
    return [ev >= 0 ? sqrt(ev) : -sqrt(-ev) for ev in evs]
end

# Input function computes the second order derivative of the RF potential
# This assumes that the good axis for radial motion does not depend on
# the ion position, otherwise we need to solve all radial modes at the same time...
function radial_modes(model::IonChainModel, r_hess)
    values = [value(ion.pos) for ion in model.ions]
    nions = length(values)
    H = zeros(nions, nions)
    for i in 1:nions
        pos = values[i]
        H[i, i] = r_hess(pos) / model.ions[i].mass
    end
    for i2 in 2:nions
        pos2 = values[i2]
        mass2 = model.ions[i2].mass
        for i1 in 1:i2 - 1
            pos1 = values[i1]
            mass1 = model.ions[i1].mass
            mass12 = sqrt(mass1 * mass2)
            term = 1 / (pos2 - pos1)^3
            H[i1, i1] -= term / mass1
            H[i2, i2] -= term / mass2
            H[i1, i2] += term / mass12
            H[i2, i1] += term / mass12
        end
    end
    evs = eigvals(H)
    return [ev >= 0 ? sqrt(ev) : -sqrt(-ev) for ev in evs]
end

function interpolate_ref_functions(fi::Ref)
    return (x->fi[](x), x->only(Interpolations.gradient(fi[], x)),
            x->only(Interpolations.hessian(fi[], x)))
end
