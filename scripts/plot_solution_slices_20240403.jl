#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
import PhoenixVoltages.Fitting
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using MAT

using NaCsPlot
using PyPlot

# Loading input data
const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end

# Utility functions
function get_rf_center(solution, xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function find_index_range(center, radius, sz)
    lb = floor(Int, center - radius)
    ub = ceil(Int, center + radius)
    lb = max(1, lb)
    ub = min(sz, ub)
    return lb:ub
end

function accumulate_electrode_slice_data!(slice_data, solution, ele, pos_index, weight)
    data = solution.data
    stride_ums_zyx = solution.stride .* 1000

    index_range_y = find_index_range(pos_index[2], 2, solution.ny)
    index_range_z = find_index_range(pos_index[3], 2, solution.nz)
    data = @view(data[index_range_z, index_range_y, :, ele])
    len = solution.nx
    fitter = Fitting.PolyFitter(4, 4)
    scales_zyx = Solutions.l_unit_um ./ stride_ums_zyx
    weight = weight / Solutions.V_unit
    for i in 1:len
        data_zy = @view(data[:, :, i])
        cache_zy = Fitting.PolyFitCache(fitter, data_zy)
        fit_zy = get(cache_zy, (pos_index[3], pos_index[2]))

        slice_data.C[i] += fit_zy[0, 0] * weight
        slice_data.Y[i] += fit_zy[0, 1] * weight * scales_zyx[2]
        slice_data.Z[i] += fit_zy[1, 0] * weight * scales_zyx[1]
        slice_data.Y2[i] += fit_zy[0, 2] * weight * scales_zyx[2]^2 * 2
        slice_data.YZ[i] += fit_zy[1, 1] * weight * scales_zyx[2] * scales_zyx[1]
        slice_data.Z2[i] += fit_zy[2, 0] * weight * scales_zyx[1]^2 * 2
    end
end

function gen_slice_data(solution, eles, voltages, pos_um)
    pos_index = get_rf_center(solution, pos_um)
    len = solution.nx
    slice_data = (C=zeros(len), Y=zeros(len), Z=zeros(len),
                  Y2=zeros(len), YZ=zeros(len), Z2=zeros(len))
    for (ele, voltage) in zip(eles, voltages)
        accumulate_electrode_slice_data!(slice_data, solution, ele,
                                         pos_index, voltage)
    end
    return slice_data
end

struct Filter
    name::String
    value::String
end

function load_mat_solution(mat, solution; filter=nothing)
    data = read(mat)
    name = get(data, "termname", nothing)
    if haskey(data, "transfer_solutions") && haskey(data, "electrode_names")
        trans_sols = data["transfer_solutions"]
        ele_names = data["electrode_names"]
        if isa(trans_sols, AbstractArray)
            if length(trans_sols) > 1
                if filter === nothing
                    error("Missing filter")
                end
                if filter.name == "#index"
                    idx = parse(Int, filter.value)
                else
                    idx = findfirst(trans_sols) do x
                        string(x[filter.name]) == filter.value
                    end
                    if idx === nothing
                        error("Solution with $(filter.name)=$(filter.value) not found")
                    end
                end
                trans_sols = trans_sols[idx]
            else
                trans_sols = trans_sols[1]
            end
        end
        pos_um = trans_sols["xpos_um"]
        eles = Int[]
        voltages = Float64[]
        for (e, v) in zip(trans_sols["electrodes"], trans_sols["voltages"])
            for ename in ele_names[e]
                push!(eles, solution.electrode_index[ename])
                push!(voltages, v)
            end
        end
        return name, pos_um, eles, voltages
    else
        error("Unknown MAT file format")
    end
end

function load_solution_file(fname, solution; filter=nothing)
    if endswith(fname, ".mat")
        return matopen(fname) do mat
            load_mat_solution(mat, solution, filter=filter)
        end
    end
    error("Unknown file format")
end

function get_filter(filter_str)
    if filter_str === nothing || isempty(filter_str)
        return nothing
    end
    filter_strs = split(filter_str, "="; limit=2)
    if length(filter_strs) < 2
        error("Invalid filter string: $(filter_str)")
    end
    return Filter(filter_strs[1], filter_strs[2])
end

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file)
plot_name, pos_um, eles, voltages =
    load_solution_file(ARGS[2], solution, filter=get_filter(get(ARGS, 3, nothing)))

slice_data = gen_slice_data(solution, eles, voltages, pos_um)

figure(figsize=[6.4 * 3, 4.8 * 2])
subplot(2, 3, 1)
plot(slice_data.C)
grid()
title("\$1\$")

subplot(2, 3, 2)
plot(slice_data.Y)
grid()
title("\$y\$")

subplot(2, 3, 3)
plot(slice_data.Z)
grid()
title("\$z\$")

subplot(2, 3, 4)
plot(slice_data.C)
grid()
title("\$y^2\$")

subplot(2, 3, 5)
plot(slice_data.Y)
grid()
title("\$yz\$")

subplot(2, 3, 6)
plot(slice_data.Z)
grid()
title("\$z^2\$")

suptitle(plot_name)

NaCsPlot.maybe_show()
