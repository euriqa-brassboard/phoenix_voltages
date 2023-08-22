#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using PhoenixVoltages
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions
using NaCsPlot
using PyPlot

const data_prefix = joinpath(@__DIR__, "../data/rf_trap_grad")
const imgs_prefix = joinpath(@__DIR__, "../imgs/rf_trap_grad")

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file)
# RF is electrode 2 (ground is 1)
const rf_data = solution.data[:, :, :, 2]

const zy_fitter = Fitting.PolyFitter(4, 4, sizes=(5, 5))

const r∇Ω_ys = Vector{Float64}(undef, solution.nx)
const r∇Ω_zs = Vector{Float64}(undef, solution.nx)

const ystride_um = solution.stride[2] * 1000
const zstride_um = solution.stride[1] * 1000

function count_idx(idxs...)
    cnt1 = 0
    cnt2 = 0
    for idx in idxs
        if idx == 1
            cnt1 += 1
        else
            cnt2 += 1
        end
    end
    return cnt1, cnt2
end

function derivative(fit, oz, oy)
    return fit[oz, oy] / zstride_um^oz / ystride_um^oy * factorial(oz) * factorial(oy)
end

for xidx in 1:solution.nx
    yidx, zidx = get(centers, xidx)
    # @show xidx, yidx, zidx
    fit_cache = Fitting.PolyFitCache(zy_fitter, @view rf_data[:, :, xidx])
    # @show (PhoenixVoltages.gradient(fit_cache, 1, zidx, yidx),
    #        PhoenixVoltages.gradient(fit_cache, 2, zidx, yidx))
    fit = get(fit_cache, (zidx, yidx))

    ∇Ω²_y = 0.0
    ∇Ω²_z = 0.0
    Ω² = 0.0

    for i in 1:2
        for j in 1:2
            cnty, cntz = count_idx(i, j)
            gij = derivative(fit, cntz, cnty)
            gijy = derivative(fit, cntz, cnty + 1)
            gijz = derivative(fit, cntz + 1, cnty)

            Ω² += gij^2
            ∇Ω²_y += 2 * gij * gijy
            ∇Ω²_z += 2 * gij * gijz
        end
    end

    r∇Ω_ys[xidx] = ∇Ω²_y / Ω² / 2
    r∇Ω_zs[xidx] = ∇Ω²_z / Ω² / 2
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "xs_um", xs_um)
    write(mat, "rel_grad", Dict("y"=>r∇Ω_ys, "z"=>r∇Ω_zs))
end

figure()
plot(xs_um, r∇Ω_ys, label="\$y\$")
plot(xs_um, r∇Ω_zs, label="\$z\$")
legend()
xlabel("X (\$\\mu m\$)")
ylabel("\$\\nabla\\Omega_{\\mathrm{RF}}/\\Omega_{\\mathrm{RF}} (\\mu m^{-1})\$")
title("\$\\Omega_{\\mathrm{RF}}\$ gradient")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)")

NaCsPlot.maybe_show()
