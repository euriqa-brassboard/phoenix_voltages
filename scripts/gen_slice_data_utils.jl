#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
import PhoenixVoltages.Fitting
using PhoenixVoltages.Potentials
using PhoenixVoltages.Mappings
using MAT

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

struct SliceData
    xpos::Vector{Float64}
    C::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    Y2::Vector{Float64}
    YZ::Vector{Float64}
    Z2::Vector{Float64}
    function SliceData(xpos)
        len = length(xpos)
        return new(xpos, zeros(len), zeros(len), zeros(len),
                   zeros(len), zeros(len), zeros(len))
    end
    Base.:*(data::SliceData, s::Number) =
        new(data.xpos, data.C .* s, data.Y .* s, data.Z .* s,
            data.Y2 .* s, data.YZ .* s, data.Z2 .* s)
    Base.:*(s::Number, data::SliceData) = data * s
    Base.:/(data::SliceData, s::Number) =
        new(data.xpos, data.C ./ s, data.Y ./ s, data.Z ./ s,
            data.Y2 ./ s, data.YZ ./ s, data.Z2 ./ s)
    Base.:\(s::Number, data::SliceData) = data / s
    Base.:+(data1::SliceData, data2::SliceData) =
        new(data1.xpos, data1.C .+ data2.C, data1.Y .+ data2.Y, data1.Z .+ data2.Z,
            data1.Y2 .+ data2.Y2, data1.YZ .+ data2.YZ, data1.Z2 .+ data2.Z2)
    Base.:-(data1::SliceData, data2::SliceData) =
        new(data1.xpos, data1.C .- data2.C, data1.Y .- data2.Y, data1.Z .- data2.Z,
            data1.Y2 .- data2.Y2, data1.YZ .- data2.YZ, data1.Z2 .- data2.Z2)
    Base.:+(data::SliceData) = data
    Base.:-(data::SliceData) =
        new(data.xpos, .-data.C, .-data.Y, .-data.Z,
            .-data.Y2, .-data.YZ, .-data.Z2)
end

function accumulate_electrode_slice_data!(slice_data, solution, ele, pos_index, weight, xidx_range)
    data = solution.data
    stride_ums_zyx = solution.stride .* 1000

    index_range_y = find_index_range(pos_index[2], 2, solution.ny)
    index_range_z = find_index_range(pos_index[3], 2, solution.nz)
    data = @view(data[index_range_z, index_range_y, xidx_range, ele])
    fitter = Fitting.PolyFitter(4, 4)
    scales_zyx = Solutions.l_unit_um ./ stride_ums_zyx
    weight = weight / Solutions.V_unit
    for i in 1:length(xidx_range)
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

function gen_slice_data(solution, eles, voltages, pos_um;
                        remove_dc=true, min_pos_um=nothing, max_pos_um=nothing)
    len = solution.nx
    idx_lo = (min_pos_um == nothing ? 1 :
        max(1, round(Int, Solutions.x_axis_to_index(solution, min_pos_um / 1000))))
    idx_hi = (max_pos_um == nothing ? len :
        min(len, round(Int, Solutions.x_axis_to_index(solution, max_pos_um / 1000))))
    pos_index = get_rf_center(solution, pos_um)
    idx_range = idx_lo:idx_hi
    idx_len = length(idx_range)
    slice_data = SliceData(Solutions.x_index_to_axis.(Ref(solution), idx_range) .* 1000)
    for (ele, voltage) in zip(eles, voltages)
        accumulate_electrode_slice_data!(slice_data, solution, ele,
                                         pos_index, voltage, idx_range)
    end
    if remove_dc
        idx_center_x = round(Int, pos_index[1])
        slice_data.C .-= slice_data.C[idx_center_x - idx_lo + 1]
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

struct Slice
    name::Union{String,Nothing}
    pos_um::Int
    eles::Vector{Int}
    voltages::Vector{Float64}
    data::SliceData
end

struct SliceLoader
    solution::Potentials.Potential
end
SliceLoader(file::AbstractString) =
    SliceLoader(Potentials.import_pillbox_64(file))

function load(loader::SliceLoader, file;
              filter=nothing, remove_dc=true, min_pos_um=nothing, max_pos_um=nothing)
    if !isa(filter, Filter)
        filter = get_filter(filter)
    end
    name, pos_um, eles, voltages =
        load_solution_file(file, loader.solution, filter=filter)
    slice_data = gen_slice_data(loader.solution, eles, voltages, pos_um;
                                remove_dc=remove_dc, min_pos_um=min_pos_um,
                                max_pos_um=max_pos_um)
    return Slice(name, pos_um, eles, voltages, slice_data)
end
