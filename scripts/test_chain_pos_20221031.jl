#!/usr/bin/julia

include("test_chain_pos_lib_20221031.jl")

using Ipopt
using PyPlot
using Test

function gen_polynomial_potential(X1, X2, X3=Ref(0.0), X4=Ref(0.0),
                                  X5=Ref(0.0), X6=Ref(0.0))
    f(x) = @evalpoly(x, 0, X1[], X2[], X3[], X4[], X5[], X6[])
    ∇f(x) = @evalpoly(x, X1[], 2 * X2[], 3 * X3[], 4 * X4[], 5 * X5[], 6 * X6[])
    ∇²f(x) = @evalpoly(x, 2 * X2[], 6 * X3[], 12 * X4[], 20 * X5[], 30 * X6[])
    return f, ∇f, ∇²f
end

function gen_model(f, ∇f, ∇²f)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    return IonChainModel(model, f, ∇f, ∇²f)
end

const X = Ref(0.0)
const X2 = Ref(2.0)
const X4 = Ref(0.0)

const chain_model = gen_model(gen_polynomial_potential(X, X2, Ref(0.0), X4)...)
const builder = ModelBuilder(chain_model)

const nions = 4

for i in 1:nions
    add_ion!(builder, i - (nions - 1) / 2, 1)
end

finalize_model!(builder)

const ion_pos = [Float64[] for i in 1:nions]
const axial_freqs = [Float64[] for i in 1:nions]
const radial_freqs = [Float64[] for i in 1:nions]

const x2s = range(0.1, 0.3, 1001)

@time for x2 in x2s
    X2[] = x2 / 2
    JuMP.optimize!(chain_model.model)
    for (ion, ions) in zip(chain_model.ions, ion_pos)
        push!(ions, value(ion.pos))
    end
    for (mode, freqs) in zip(axial_modes(chain_model), axial_freqs)
        push!(freqs, mode)
    end
    for (mode, freqs) in zip(radial_modes(chain_model, x->4), radial_freqs)
        push!(freqs, mode)
    end
    update_init_pos!(chain_model)
end

figure()
for ions in ion_pos
    plot(x2s, ions)
end
grid()

figure()
for freqs in axial_freqs
    plot(x2s, freqs)
end
grid()

figure()
for freqs in radial_freqs
    plot(x2s, freqs)
end
grid()

show()
