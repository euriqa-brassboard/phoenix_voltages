#!/usr/bin/julia

include("test_chain_pos_lib_20221031.jl")

using NLopt
using Ipopt
using PyPlot

# const model = Model(NLopt.Optimizer)
# set_optimizer_attribute(model, "algorithm", :LD_MMA)
const model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

const X = Ref(0.0)
const X2 = Ref(0.0)
const X4 = Ref(0.2)

f(x) = (X2[] + X4[] * x^2) * x^2 + X[] * x
∇f(x) = (2 * X2[] + 4 * X4[] * x^2) * x + X[]
∇²f(x) = 2 * X2[] + 12 * X4[] * x^2

const chain_model = IonChainModel(model, f, ∇f, ∇²f)
const builder = ModelBuilder(chain_model)

const nions = 30

for i in 1:nions
    add_ion!(builder, i - (nions - 1) / 2, 1)
end

finalize_model!(builder)

const ion_pos = [Float64[] for i in 1:nions]
const axial_freqs = [Float64[] for i in 1:nions]
const radial_freqs = [Float64[] for i in 1:nions]

const x2s = range(0, -10, 1001)

@time for x2 in x2s
    X2[] = x2
    JuMP.optimize!(model)
    for (ion, ions) in zip(chain_model.ions, ion_pos)
        push!(ions, value(ion.pos))
    end
    for (mode, freqs) in zip(axial_modes(chain_model), axial_freqs)
        push!(freqs, mode)
    end
    for (mode, freqs) in zip(radial_modes(chain_model, x->5000), radial_freqs)
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
