#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.ProcessSolution
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using NaCsPlot
using PyPlot
using MAT
using LsqFit
using LinearAlgebra

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return ProcessSolution.CenterTracker(read(mat, "zy_index"))
end
const short_map = Dict{String,String}()

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = ProcessSolution.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)

function get_rf_center(xpos_um)
    xidx = ProcessSolution.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function get_compensate_terms1(xpos_um)
    @show xpos_um
    return ProcessSolution.get_compensate_terms1(fits_cache, get_rf_center(xpos_um))
end

function get_compensate_fits1(xpos_um, electrode)
    @show xpos_um
    pos = get_rf_center(xpos_um)
    return get(fits_cache, electrode, (pos[3], pos[2], pos[1]))
end

function get_all_names(ids)
    names = String[]
    for id in ids
        append!(names, solution.electrode_names[id])
    end
    sort!(names)
end

# const xpos_ums = -3220:1330
# const xpos_ums = -1800:-1600
const xpos_um1 = -1712
const xpos_um2 = -1702

const electrodes = ["S" .* string.(0:11); "O0"; "O1";
                    "Q44"; "Q45"; "Q0"; "Q1"; "Q2"; "Q3"]
const electrode_ids = sort!(unique(getindex.(Ref(solution.electrode_index),
                                             electrodes)))
const comp_fits1 = [get_compensate_fits1(xpos_um1, ele) for ele in electrode_ids]
const comp_fits2 = [get_compensate_fits1(xpos_um2, ele) for ele in electrode_ids]

const comp_terms1 = [ProcessSolution.get_compensate_terms1(fit, solution.stride .* 1000) for fit in comp_fits1]
const comp_terms2 = [ProcessSolution.get_compensate_terms1(fit, solution.stride .* 1000) for fit in comp_fits2]

function interpolate_fit(fit1, fit2, x)
    return (dx=fit1.dx .* (1 - x) .+ fit2.dx .* x,
            dy=fit1.dy, # .* (1 - x) .+ fit2.dy .* x,
            dz=fit1.dz, # .* (1 - x) .+ fit2.dz .* x,

            xy=fit1.xy, # .* (1 - x) .+ fit2.xy .* x,
            yz=fit1.yz, # .* (1 - x) .+ fit2.yz .* x,
            zx=fit1.zx .* (1 - x) .+ fit2.zx .* x,

            z2=fit1.z2, # .* (1 - x) .+ fit2.z2 .* x,
            x2=fit1.x2, # .* (1 - x) .+ fit2.x2 .* x,
            x3=fit1.x3, # .* (1 - x) .+ fit2.x3 .* x,
            x4=fit1.x4, # .* (1 - x) .+ fit2.x4 .* x,
            )
end

# function interpolate_fit(fit1, fit2, x)
#     return (dx=fit1.dx, # .* (1 - x) .+ fit2.dx .* x,
#             dy=fit1.dy .* (1 - x) .+ fit2.dy .* x,
#             dz=fit1.dz .* (1 - x) .+ fit2.dz .* x,

#             xy=fit1.xy .* (1 - x) .+ fit2.xy .* x,
#             yz=fit1.yz .* (1 - x) .+ fit2.yz .* x,
#             zx=fit1.zx, # .* (1 - x) .+ fit2.zx .* x,

#             z2=fit1.z2 .* (1 - x) .+ fit2.z2 .* x,
#             x2=fit1.x2 .* (1 - x) .+ fit2.x2 .* x,
#             x3=fit1.x3 .* (1 - x) .+ fit2.x3 .* x,
#             x4=fit1.x4 .* (1 - x) .+ fit2.x4 .* x,
#             )
# end

function interpolate_fits(fits1, fits2, x)
    return [interpolate_fit(fit1, fit2, x) for (fit1, fit2) in zip(fits1, fits2)]
end

interpolate_fits(x) = interpolate_fits(comp_terms1, comp_terms2, x)

function solve_terms2(terms)
    nfits = length(terms)
    coefficient = Matrix{Float64}(undef, 10, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(terms[i])
    end
    X = coefficient \ Matrix(I, 10, 10)
    @assert size(X, 2) == 10
    return (dx=X[:, 1], dy=X[:, 2], dz=X[:, 3],
            xy=X[:, 4], yz=X[:, 5], zx=X[:, 6],
            z2=X[:, 7], x2=X[:, 8], x3=X[:, 9], x4=X[:, 10])
end
solve_interpolate(x) = solve_terms2(interpolate_fits(x))

const comp_terms_sol1 = solve_terms2(comp_terms1)
const comp_terms_sol2 = solve_terms2(comp_terms2)


# @show comp_terms_sol1
# @show comp_terms_sol2

@show maximum(abs.(comp_terms_sol1.x2))
@show maximum(abs.(comp_terms_sol2.x2))

@show maximum(abs.(solve_interpolate(0).x2))
@show maximum(abs.(solve_interpolate(0.5).x2))
@show maximum(abs.(solve_interpolate(1).x2))

const interpolate_xs = range(-0.1, 1.1, 1000)

figure()
plot(interpolate_xs, [maximum(abs.(solve_interpolate(x).x2)) for x in interpolate_xs])
grid()

# const comp_terms = [get_compensate_terms1(xpos_um) for xpos_um in xpos_ums]
# for term in comp_terms
#     println(join(get_all_names(term[1]), "."))
# end
# const comp_files = [ProcessSolution.compensation_to_file(solution, mapfile, term...)
#                     for term in comp_terms]

NaCsPlot.maybe_show()
