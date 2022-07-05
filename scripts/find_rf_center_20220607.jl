#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using NaCsPlot
using PyPlot
using MAT

const data_prefix = joinpath(@__DIR__, "../data/rf_center")
const imgs_prefix = joinpath(@__DIR__, "../imgs/rf_center")

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file)
# RF is electrode 2 (ground is 1)
const centers = Solutions.find_all_flat_points(solution.data[:, :, :, 2])
const centers_um = similar(centers)
for i in 1:solution.nx
    centers_um[i, 1] = z_index_to_axis(solution, centers[i, 1]) * 1000
    centers_um[i, 2] = y_index_to_axis(solution, centers[i, 2]) * 1000
end

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "zy_index", centers)
    write(mat, "zy_um", centers_um)
end

const residual_grad = similar(centers)
# This is mainly a test for the shift function...
for i in 1:solution.nx
    data = solution.data[:, :, i, 2]
    fitter = Fitting.PolyFitter(3, 3)
    cache = Fitting.PolyFitCache(fitter, data)
    fit = get(cache, (centers[i, 1], centers[i, 2]))
    residual_grad[i, 1] = fit[1, 0]
    residual_grad[i, 2] = fit[0, 1]
    # residual_grad[i, 1] = PhoenixVoltages.gradient(cache, 1, centers[i, 1], centers[i, 2])
    # residual_grad[i, 2] = PhoenixVoltages.gradient(cache, 2, centers[i, 1], centers[i, 2])
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

figure()
plot(xs_um, residual_grad[:, 1], label="Z")
plot(xs_um, residual_grad[:, 2], label="Y")
xlabel("X (\$\\mu\$m)")
title("Residual gradient")
grid()

figure(figsize=[6.4 * 2, 4.8])
subplot(1, 2, 1)
plot(xs_um, centers_um[:, 1])
xlabel("X (\$\\mu\$m)")
ylabel("Z (\$\\mu\$m)")
title("RF null Z position")
grid()

subplot(1, 2, 2)
plot(xs_um, centers_um[:, 2])
xlabel("X (\$\\mu\$m)")
ylabel("Y (\$\\mu\$m)")
title("RF null Y position")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)_zy")

NaCsPlot.maybe_show()
