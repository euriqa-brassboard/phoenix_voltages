#!/usr/bin/julia

include("test_chain_pos_lib_20221031.jl")

using Ipopt
using PyPlot

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

const nions = 33

for i in 1:nions
    add_ion!(builder, i - (nions - 1) / 2, 1)
end

finalize_model!(builder)

const ion_pos = [Float64[] for i in 1:nions]
const axial_freqs = [Float64[] for i in 1:nions]
const radial_freqs = [Float64[] for i in 1:nions]

function interpolate(x, xs, ys)
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end
    idx = searchsortedfirst(xs, x)
    x0 = xs[idx - 1]
    x1 = xs[idx]
    y0 = ys[idx - 1]
    y1 = ys[idx]
    x = (x - x0) / (x1 - x0)
    return y1 * x + y0 * (1 - x)
end

const X1_0 = 0.094681
const X2_0 = 0.02
const X4_0 = 0.015e-3

function calc_x2(a)
    return a * X2_0
end

function calc_x4(a)
    return interpolate(a, [-0.3, 0, 1], (1.0, 1.0, 0.0)) * X4_0
end
function calc_x1(a)
    return interpolate(a, [-0.3, 0], (1.0, 0.0)) * X1_0
end

const as = range(1, -0.5, 1001)

@time for a in as
    X[] = calc_x1(a)
    X2[] = calc_x2(a) / 2
    X4[] = calc_x4(a) / 2 / 3 / 4
    JuMP.optimize!(chain_model.model)
    for (ion, ions) in zip(chain_model.ions, ion_pos)
        push!(ions, value(ion.pos))
    end
    for (mode, freqs) in zip(axial_modes(chain_model), axial_freqs)
        push!(freqs, mode)
    end
    for (mode, freqs) in zip(radial_modes(chain_model, x->6)[1], radial_freqs)
        push!(freqs, mode)
    end
    update_init_pos!(chain_model)
end

figure()
for ions in ion_pos
    plot(as, ions)
end
grid()

figure()
for freqs in axial_freqs
    plot(as, freqs)
end
grid()

figure()
for freqs in radial_freqs
    plot(as, freqs)
end
grid()

show()

# X1 range for merging a single ion
# 32: 0.094680 [0.094681, 0.106609] 0.106610
# 16: 0.074121 [0.074122, 0.090556] 0.090557
#  8: 0.052261 [0.052262, 0.076170] 0.076171
#  2: -.048342 [-.048341, 0.048341] 0.048342

# Ion number range for each X1
# 0.048341: [2, 7]
# 0.052262: [3, 8]
# 0.074122: [8, 16]
# 0.076170: [8, 17]
# 0.090556: [16, 27]
# 0.094681: [20, 32]
