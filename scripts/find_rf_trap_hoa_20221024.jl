#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using PhoenixVoltages
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions
using NaCsPlot
using PyPlot

const data_prefix = joinpath(@__DIR__, "../data/rf_trap_hoa")
const imgs_prefix = joinpath(@__DIR__, "../imgs/rf_trap_hoa")

const centers = matopen(joinpath(@__DIR__, "../data/rf_center_hoa.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, is_hoa=true)
# RF is electrode 2 (ground is 1)
const rf_data = solution.data[:, :, :, 2]

const zy_fitter = Fitting.PolyFitter(4, 4, sizes=(5, 5))

const rf_y2s = Vector{Float64}(undef, solution.nx)
const rf_yzs = Vector{Float64}(undef, solution.nx)
const rf_z2s = Vector{Float64}(undef, solution.nx)

const trap_y2s = Vector{Float64}(undef, solution.nx)
const trap_yzs = Vector{Float64}(undef, solution.nx)
const trap_z2s = Vector{Float64}(undef, solution.nx)

const trap_min = Vector{Float64}(undef, solution.nx)
const trap_max = Vector{Float64}(undef, solution.nx)

const ystride_m = solution.stride[2] / 1000
const zstride_m = solution.stride[1] / 1000

const e = 1.60217663e-19 # C
const Ω_trap = 2π * 46.25e6 # Hz
const NA = 6.02214076e23
const m_Yb171 = 170.9363315e-3 / NA # kg
# the last 2pi converts angular frequency back to real frequency
const scale_rf = 2 * (e / m_Yb171 / Ω_trap / 2π)^2

for xidx in 1:solution.nx
    yidx, zidx = get(centers, xidx)
    # @show xidx, yidx, zidx
    fit_cache = Fitting.PolyFitCache(zy_fitter, @view rf_data[:, :, xidx])
    # @show (PhoenixVoltages.gradient(fit_cache, 1, zidx, yidx),
    #        PhoenixVoltages.gradient(fit_cache, 2, zidx, yidx))
    fit = get(fit_cache, (zidx, yidx))
    y2 = fit[0, 2] / ystride_m^2
    yz = fit[1, 1] / zstride_m / ystride_m / 2
    z2 = fit[2, 0] / zstride_m^2
    rf_y2s[xidx] = y2 # 1/m^2
    rf_yzs[xidx] = yz # 1/m^2
    rf_z2s[xidx] = z2 # 1/m^2

    # Compute the matrix the determines the trap potential shape
    α2 = [y2 yz; yz z2]^2 .* scale_rf # Hz^2 / V^2
    a = α2[1, 1]
    b = (α2[2, 1] + α2[1, 2]) / 2
    c = α2[2, 2]

    trap_y2s[xidx] = a
    trap_yzs[xidx] = b
    trap_z2s[xidx] = c

    Δ = sqrt((a - c)^2 + 4 * b^2)

    trap_min[xidx] = sqrt((a + c - Δ) / 2) # Hz/V
    trap_max[xidx] = sqrt((a + c + Δ) / 2) # Hz/V
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "xs_um", xs_um)
    write(mat, "field", Dict("y2"=>rf_y2s, "yz"=>rf_yzs, "z2"=>rf_z2s))
    write(mat, "trap", Dict("y2"=>trap_y2s, "yz"=>trap_yzs, "z2"=>trap_z2s))
    write(mat, "frequency", Dict("min"=>trap_min, "max"=>trap_max))
end

figure()
plot(xs_um, rf_y2s ./ 1e6, label="\$y^2\$")
plot(xs_um, rf_yzs ./ 1e6, label="\$2yz\$")
plot(xs_um, rf_z2s ./ 1e6, label="\$z^2\$")
legend()
xlabel("X (\$\\mu m\$)")
ylabel("\$mm^{-2}\$")
title("RF Field Curvature")
ylim([-25, 25])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_field")

figure()
plot(xs_um, rf_y2s ./ 1e6, label="\$y^2\$")
plot(xs_um, .-rf_z2s ./ 1e6, label="\$-z^2\$")
legend()
xlabel("X (\$\\mu m\$)")
ylabel("\$mm^{-2}\$")
title("RF Field Curvature")
grid()
ylim([22.5, 24.5])
NaCsPlot.maybe_save("$(imgs_prefix)_field_zhj")

figure()
plot(xs_um, .-(rf_y2s .+ rf_z2s) ./ 1e6, label="\$x^2\$")
legend()
xlabel("X (\$\\mu m\$)")
ylabel("\$mm^{-2}\$")
title("RF Field Axial Curvature")
ylim([-0.2, 0.2])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_field_axial")

figure()
plot(xs_um, trap_y2s ./ 1e6, label="\$y^2\$")
plot(xs_um, trap_yzs ./ 1e6, label="\$2yz\$")
plot(xs_um, trap_z2s ./ 1e6, label="\$z^2\$")
legend()
xlabel("X (\$\\mu m\$)")
ylabel("\$kHz^2/V^2\$")
title("Ion Trap Potential Curvature")
ylim([-3, 120])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_potential")

figure()
plot(xs_um, trap_min ./ 1e3, label="minor")
plot(xs_um, trap_max ./ 1e3, label="major")
legend()
xlabel("X (\$\\mu m\$)")
ylabel("\$kHz/V\$")
title("Secular Frequency")
ylim([9.8, 10.7])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_frequency")

NaCsPlot.maybe_show()
