#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using PhoenixVoltages
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions
using NaCsPlot
using PyPlot

const data_prefix = joinpath(@__DIR__, "../data/rf_axial_peregrine")
const imgs_prefix = joinpath(@__DIR__, "../imgs/rf_axial_peregrine")

const centers = matopen(joinpath(@__DIR__, "../data/rf_center_peregrine.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file)
# RF is electrode 2 (ground is 1)
const rf_data = solution.data[:, :, :, 2]

const fitter = Fitting.PolyFitter(3, 3, 4, sizes=(5, 5, 15))
const fit_cache = Fitting.PolyFitCache(fitter, rf_data)

const trap_axial = Vector{Float64}(undef, solution.nx)

const trap_x2 = Vector{Float64}(undef, solution.nx)
const trap_y2 = Vector{Float64}(undef, solution.nx)
const trap_z2 = Vector{Float64}(undef, solution.nx)

const xstride_m = solution.stride[3] / 1000

const xstride_mm = solution.stride[3]
const ystride_mm = solution.stride[2]
const zstride_mm = solution.stride[1]

for xidx in 1:solution.nx
    yidx, zidx = get(centers, xidx)
    # @show xidx, yidx, zidx
    fit = get(fit_cache, (zidx, yidx, Float64(xidx)))
    # fit[0, 0, 1] is in 1/um
    trap_axial[xidx] = fit[0, 0, 1] / xstride_m

    trap_x2[xidx] = fit[0, 0, 2] / xstride_mm^2
    trap_y2[xidx] = fit[0, 2, 0] / ystride_mm^2
    trap_z2[xidx] = fit[2, 0, 0] / zstride_mm^2
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "xs_um", xs_um)
    write(mat, "field", trap_axial)
end

figure()
plot(xs_um, trap_axial)
xlabel("X (\$\\mu m\$)")
ylabel("Axial RF field (\$m^{-1}\$)")
ylim([-7, 7])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)")

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
ylim([-0.1, 0.1])
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_0")

NaCsPlot.maybe_show()
