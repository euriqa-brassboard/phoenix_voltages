#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using PhoenixVoltages
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions
using NaCsPlot
using PyPlot

const data_prefix = joinpath(@__DIR__, "../data/rf_axial")
const imgs_prefix = joinpath(@__DIR__, "../imgs/rf_axial")

const centers = Solutions.CenterTracker()

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file)
# RF is electrode 2 (ground is 1)
const rf_data = solution.data[:, :, :, 2]

const fitter = Fitting.PolyFitter(4, 4, 4, sizes=(5, 5, 15))
const fit_cache = Fitting.PolyFitCache(fitter, rf_data)

const trap_offset = Vector{Float64}(undef, solution.nx)
const trap_offset_min = Vector{Float64}(undef, solution.nx)
const trap_offset_max = Vector{Float64}(undef, solution.nx)
const trap_axial = Vector{Float64}(undef, solution.nx)

const trap_x2 = Vector{Float64}(undef, solution.nx)
const trap_y2 = Vector{Float64}(undef, solution.nx)
const trap_z2 = Vector{Float64}(undef, solution.nx)

const trap_offset_pseudo = Vector{Float64}(undef, solution.nx)
const trap_axial_pseudo = Vector{Float64}(undef, solution.nx)

const xstride_m = solution.stride[3] / 1000

const xstride_mm = solution.stride[3]
const ystride_mm = solution.stride[2]
const zstride_mm = solution.stride[1]

const e = 1.60217663e-19 # C
const Ω_trap = 2π * 46.25e6 # Hz
const NA = 6.02214076e23
const m_Yb171 = 170.9363315e-3 / NA # kg
const scale_rf_pseudo = e / 4 / m_Yb171 / Ω_trap^2 # V / (V/m)^2

for xidx in 1:solution.nx
    yidx, zidx = get(centers, xidx)
    # @show xidx, yidx, zidx
    fit = get(fit_cache, (zidx, yidx, Float64(xidx)))
    trap_offset[xidx] = fit[0, 0, 0]
    trap_offset_min[xidx] = minimum(rf_data[:, :, xidx])
    trap_offset_max[xidx] = maximum(rf_data[:, :, xidx])
    field = fit[0, 0, 1] / xstride_m # m^-1
    curvature = fit[0, 0, 2] / xstride_m^2 # m^-2
    trap_axial[xidx] = field

    trap_offset_pseudo[xidx] = scale_rf_pseudo * field^2 # V / V^2
    trap_axial_pseudo[xidx] = scale_rf_pseudo * 4 * field * curvature # V/m / V^2

    trap_x2[xidx] = fit[0, 0, 2] / xstride_mm^2
    trap_y2[xidx] = fit[0, 2, 0] / ystride_mm^2
    trap_z2[xidx] = fit[2, 0, 0] / zstride_mm^2
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "xs_um", xs_um)
    write(mat, "offset", trap_offset)
    write(mat, "field", trap_axial)
    write(mat, "pseudo_offset", trap_offset_pseudo)
    write(mat, "pseudo_field", trap_axial_pseudo)
    write(mat, "x2", trap_x2)
    write(mat, "y2", trap_y2)
    write(mat, "z2", trap_z2)
end

figure()
plot(xs_um, trap_axial)
xlabel("X (\$\\mu m\$)")
ylabel("Axial RF field (\$m^{-1}\$)")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)")

figure()
plot(xs_um, trap_axial)
xlabel("X (\$\\mu m\$)")
ylabel("Axial RF field (\$m^{-1}\$)")
ylim([-7, 7])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_zoomin")

figure()
plot(xs_um, trap_offset, "C1-")
plot(xs_um, trap_offset_min, "C0--")
plot(xs_um, trap_offset_max, "C2--")
xlabel("X (\$\\mu m\$)")
ylabel("RF potential")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_offset")

figure()
plot(xs_um, trap_x2)
plot(xs_um, trap_y2)
plot(xs_um, trap_z2)
xlabel("X (\$\\mu m\$)")
ylabel("RF field curvature (\$mm^{-2}\$)")
ylim([-25, 25])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_xyz2")

figure()
plot(xs_um, trap_x2 .+ trap_y2 .+ trap_z2)
xlabel("X (\$\\mu m\$)")
ylabel("Violation (\$mm^{-2}\$)")
ylim([-0.01, 0.01])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_0")

figure()
plot(xs_um, trap_offset_pseudo)
xlabel("X (\$\\mu m\$)")
ylabel("Axial pseudopotential (\$V / V^2\$)")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_pseudo")

figure()
plot(xs_um, trap_axial_pseudo)
xlabel("X (\$\\mu m\$)")
ylabel("Axial pseudopotential gradient (\$V/m / V^2\$)")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_pseudo_grad")

NaCsPlot.maybe_show()
