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

const fitter = Fitting.PolyFitter(3, 3, 4, sizes=(5, 5, 15))
const fit_cache = Fitting.PolyFitCache(fitter, rf_data)

const trap_axial = Vector{Float64}(undef, solution.nx)

const xstride_m = solution.stride[3] / 1000

for xidx in 1:solution.nx
    yidx, zidx = get(centers, xidx)
    # @show xidx, yidx, zidx
    fit = get(fit_cache, (zidx, yidx, Float64(xidx)))
    # fit[0, 0, 1] is in 1/um
    trap_axial[xidx] = fit[0, 0, 1] / xstride_m
end

const xs_um = x_index_to_axis.(Ref(solution), 1:solution.nx) .* 1000

matopen("$(data_prefix).mat", "w") do mat
    write(mat, "xs_um", xs_um)
    write(mat, "field", trap_axial)
end

figure()
plot(xs_um, trap_axial)
xlabel("X (\$\\mu m\$)")
ylabel("Axial RF field \$m^{-1}\$")
grid()
NaCsPlot.maybe_save("$(imgs_prefix)")

NaCsPlot.maybe_show()
