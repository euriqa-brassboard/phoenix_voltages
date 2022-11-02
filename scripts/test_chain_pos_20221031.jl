#!/usr/bin/julia

include("test_chain_pos_lib_20221031.jl")

using NLopt
using Ipopt
using PyPlot

# const model = Model(NLopt.Optimizer)
# set_optimizer_attribute(model, "algorithm", :LD_MMA)
const model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

const X2 = Ref(0.0)
const X4 = Ref(0.2)

f(x) = (X2[] + X4[] * x^2) * x^2 + x * 2
∇f(x) = (2 * X2[] + 4 * X4[] * x^2) * x + 2
∇²f(x) = 2 * X2[] + 12 * X4[] * x^2

const chain_model = IonChainModel(model, f, ∇f, ∇²f)
const builder = ModelBuilder(chain_model)

const nions = 31

for i in 1:nions
    add_ion!(builder, i - (nions - 1) / 2, 1)
end

finalize_model!(builder)

const ion_pos = [Float64[] for i in 1:nions]

const x2s = range(0, -10, 1001)

@time for x2 in x2s
    X2[] = x2
    JuMP.optimize!(model)
    for ((ion, charge), ions) in zip(chain_model.ions, ion_pos)
        push!(ions, value(ion))
    end
    update_init_pos!(chain_model)
end

figure()
for ions in ion_pos
    plot(x2s, ions)
end
grid()

show()
