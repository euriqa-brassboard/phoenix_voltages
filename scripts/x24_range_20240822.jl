#!/usr/bin/julia

include("test_chain_pos_lib_20221031.jl")

using Ipopt
using PyPlot
using NaCsPlot

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
const X2 = Ref(0.0)
const X4 = Ref(0.0)

const chain_model = gen_model(gen_polynomial_potential(X, X2, Ref(0.0), X4)...)
const builder = ModelBuilder(chain_model)

const nions = 32

for i in 1:nions
    add_ion!(builder, i - (nions - 1) / 2, 1)
end

finalize_model!(builder)

const x2s = range(0, 0.036, 201)
const x4s = range(0, 0.001, 201) .* 1000

const min_radials = zeros(length(x2s), length(x4s))
const max_distances = zeros(length(x2s), length(x4s))

const prefix = joinpath(@__DIR__, "../imgs/x24_range_20240822")

@time for (i2, x2) in enumerate(x2s)
    X2[] = x2 / 2
    for (i4, x4) in enumerate(x4s)
        X4[] = x4 / 2 / 3 / 4 / 1000
        JuMP.optimize!(chain_model.model)
        minpos, maxpos = extrema(value(ion.pos) for ion in chain_model.ions)
        max_distances[i2, i4] = maxpos - minpos
        min_radials[i2, i4] = minimum(radial_modes(chain_model, x->6)[1])
        update_init_pos!(chain_model)
    end
end

figure()
imshow(min_radials, origin="lower", aspect="auto", vmin=0,
       extent=(first(x4s) - step(x4s) / 2, last(x4s) + step(x4s) / 2,
               first(x2s) - step(x2s) / 2, last(x2s) + step(x2s) / 2))
xlim([first(x4s), last(x4s)])
ylim([first(x2s), last(x2s)])
xlabel("x4 * 1000")
ylabel("x2")
colorbar()
grid()
title("Minimum Radial Frequency (MHz)")
NaCsPlot.maybe_save("$(prefix)_radial_freq")

figure()
imshow(max_distances, origin="lower", aspect="auto", vmax=60, vmin=40,
       extent=(first(x4s) - step(x4s) / 2, last(x4s) + step(x4s) / 2,
               first(x2s) - step(x2s) / 2, last(x2s) + step(x2s) / 2))
xlim([first(x4s), last(x4s)])
ylim([first(x2s), last(x2s)])
xlabel("x4 * 1000")
ylabel("x2")
colorbar()
grid()
title("Chain Length (\$\\mu m\$)")
NaCsPlot.maybe_save("$(prefix)_chain_length")

NaCsPlot.maybe_show()
