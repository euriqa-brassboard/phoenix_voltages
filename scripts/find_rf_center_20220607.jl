#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.VoltageSolutions
using NaCsPlot
using PyPlot
using MAT

const data_prefix = joinpath(@__DIR__, "../data/rf_center")
const imgs_prefix = joinpath(@__DIR__, "../imgs/rf_center")

const solution_file = ARGS[1]
const solution = VoltageSolutions.import_pillbox_64(solution_file)
# RF is electrode 2 (ground is 1)
const centers = ProcessSolution.find_all_flat_points(solution.data[:, :, :, 2])
const centers_um = similar(centers)
for i in 1:solution.nx
    centers_um[i, 1] = z_index_to_axis(solution, centers[i, 1]) * 1000
    centers_um[i, 2] = y_index_to_axis(solution, centers[i, 2]) * 1000
end

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "zy_index", centers)
    write(mat, "zy_um", centers_um)
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

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
