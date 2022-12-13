#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
import PhoenixVoltages.Mappings
using PhoenixVoltages.Potentials
using MAT
using LinearAlgebra

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_red_202212.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = Solutions.compensate_fitter3(solution)

const prefix = joinpath(@__DIR__, "../data/compensate_red2_20221213")

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

# Copy of Solutions.get_compensate_coeff1 but using fewer electrodes
function get_compensate_coeff1(cache::Potentials.FitCache, pos::NTuple{3})
    # pos is in xyz index

    x_coord = Solutions.x_index_to_axis(cache.solution, pos[1]) .* 1000
    ele_select = Mappings.find_electrodes(cache.solution.electrode_index,
                                          x_coord, min_num=12, min_dist=210)
    ele_select = sort!(collect(ele_select))
    fits = [get(cache, e, (pos[3], pos[2], pos[1])) for e in ele_select]

    # Change stride to um in unit
    stride_um = cache.solution.stride .* 1000
    nfits = length(fits)
    coefficient = Matrix{Float64}(undef, 10, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(Solutions.get_compensate_terms1(fits[i], stride_um))
    end
    return ele_select, coefficient
end

function get_compensate_coeff_data(xpos_um)
    @show xpos_um
    eles, coeff = get_compensate_coeff1(fits_cache, get_rf_center(xpos_um))
    x0 = coeff \ Matrix(I, 10, 10)
    B = qr(coeff').Q[:, (10 + 1):end]
    return Dict("electrodes"=>eles,
                "coefficients"=>coeff,
                "solution"=>x0,
                "free_solution"=>B,
                "xpos_um"=>xpos_um)
end

const xpos_ums = -500:500
const coeff_data = [@time(get_compensate_coeff_data(xpos_um)) for xpos_um in xpos_ums]

matopen("$(prefix).mat", "w") do mat
    write(mat, "data", coeff_data)
    write(mat, "electrode_names", solution.electrode_names)
end
