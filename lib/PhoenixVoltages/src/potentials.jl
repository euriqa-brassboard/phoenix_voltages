#!/usr/bin/julia

module Potentials

import ..Fitting, ..get_single

struct RawPotential
    electrodes::Int
    nx::Int
    ny::Int
    nz::Int
    vsets::Int
    xaxis::NTuple{3,Float64}
    yaxis::NTuple{3,Float64}
    stride::NTuple{3,Float64}
    origin::NTuple{3,Float64}
    electrodemapping::Vector{Int}
    data::Array{Float64,4}
end

"""
    import_pillbox_v0(filename) -> (header, data)

Imports voltage array files of format V0, here `filename` is the name of the file
to be read. Voltage arrays contain potential data for one electrode at 1V
and all other electrodes at 0V on a 3D grid of points.

The data is returned in the header and data fields.
header contains the fields 'electrodes' for the number of potentials
for different electrodes,
'nx', 'ny', 'nz' are the number of samples in x-, y-, and z- direction.
'origin' is one end point of the 3D grid,
'stride' is the stepsize in the 3 dimensions.
"""
function import_pillbox_v0_raw(filename)
    open(filename) do fh
        read(fh, Int32) # discard
        electrodes = Int(read(fh, Int32))
        nx = Int(read(fh, Int32))
        ny = Int(read(fh, Int32))
        nz = Int(read(fh, Int32))
        vsets = Int(read(fh, Int32))
        stride = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
        origin = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
        # I have no idea what's stored in these
        read(fh, Int32)
        read(fh, Int32)
        electrodemapping = Vector{Int}(undef, electrodes)
        for i in 1:electrodes
            electrodemapping[i] = read(fh, Int32)
        end
        databytes = read(fh)
        if length(databytes) != electrodes * nx * ny * nz * sizeof(Float64)
            error("Did not find the right number of samples")
        end
        data = Array{Float64}(undef, nz, ny, nx, electrodes)
        copyto!(data, reinterpret(Float64, databytes))
        return RawPotential(electrodes, nx, ny, nz, vsets,
                            (1000, 0, 0), (0, 1000, 0),
                            stride, origin, electrodemapping, data)
    end
end

"""
    import_pillbox_v1(filename) -> (header, data)

Imports voltage array files of format V1, the current format,
here filename is the name of the file to be read.
Voltage arrays contain potential data for one electrode at 1V
and all other electrodes at 0V on a 3D grid of points.

The data is returned in the header and data fields.
header contains the fields 'electrodes' for the number of potentials
for different electrodes, 'nx', 'ny', 'nz' are the number of samples
in x-, y-, and z- direction.
'origin' is one end point of the 3D grid,
'stride' is the stepsize in the 3 dimensions.
"""
function import_pillbox_v1_raw(filename)
    open(filename) do fh
        read(fh, Int32) # discard
        electrodes = Int(read(fh, Int32))
        nx = Int(read(fh, Int32))
        ny = Int(read(fh, Int32))
        nz = Int(read(fh, Int32))
        vsets = Int(read(fh, Int32))
        xaxis = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64))
        yaxis = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64))
        stride = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
        origin = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
        # I have no idea what's stored in these
        read(fh, Int32)
        read(fh, Int32)
        electrodemapping = Vector{Int}(undef, electrodes)
        for i in 1:electrodes
            electrodemapping[i] = read(fh, Int32)
        end
        databytes = read(fh)
        if length(databytes) != electrodes * nx * ny * nz * sizeof(Float64)
            error("Did not find the right number of samples")
        end
        data = Array{Float64}(undef, nz, ny, nx, electrodes)
        copyto!(data, reinterpret(Float64, databytes))
        return RawPotential(electrodes, nx, ny, nz, vsets, xaxis, yaxis,
                            stride, origin, electrodemapping, data)
    end
end

"""
    import_pillbox_64(filename) -> (header, data)

Imports voltage array files in 64 bit format,
here filename is the name of the file to be read.
Voltage arrays contain potential data for one electrode at 1V
and all other electrodes at 0V on a 3D grid of points.

The data is returned in the header and data fields.
header contains the fields 'electrodes' for the number of potentials
for different electrodes, 'nx', 'ny', 'nz' are the number of samples
in x-, y-, and z- direction.
'origin' is one end point of the 3D grid,
'stride' is the stepsize in the 3 dimensions.
"""
function import_pillbox_64_raw(filename)
    open(filename) do fh
        read(fh, Int64) # discard
        electrodes = Int(read(fh, Int64))
        nx = Int(read(fh, Int64))
        ny = Int(read(fh, Int64))
        nz = Int(read(fh, Int64))
        vsets = Int(read(fh, Int64))
        xaxis = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64))
        yaxis = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64))
        stride = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
        origin = 1000 .* (read(fh, Float64), read(fh, Float64), read(fh, Float64)) # Use mm instead of m
        # I have no idea what's stored in these
        read(fh, Int64)
        read(fh, Int64)
        electrodemapping = Vector{Int}(undef, electrodes)
        for i in 1:electrodes
            electrodemapping[i] = read(fh, Int64)
        end
        databytes = read(fh)
        if length(databytes) != electrodes * nx * ny * nz * sizeof(Float64)
            error("Did not find the right number of samples")
        end
        data = Array{Float64}(undef, nz, ny, nx, electrodes)
        copyto!(data, reinterpret(Float64, databytes))
        return RawPotential(electrodes, nx, ny, nz, vsets, xaxis, yaxis,
                            stride, origin, electrodemapping, data)
    end
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        export $(Symbol(name, "_index_to_axis"))
        $(Symbol(name, "_index_to_axis"))(sol::RawPotential, i) = (i - 1) * sol.stride[$i] + sol.origin[$i]
        export $(Symbol(name, "_axis_to_index"))
        $(Symbol(name, "_axis_to_index"))(sol::RawPotential, a) = (a - sol.origin[$i]) / sol.stride[$i] + 1
    end
end

const raw_electrode_names = ["GND"; "RF";
                             "L" .* string.(0:9);
                             "O" .* string.(0:1);
                             "Q" .* string.(0:65);
                             "S" .* string.(0:11);]
const raw_electrode_index = Dict{String,Int}()
for i in 1:length(raw_electrode_names)
    raw_electrode_index[raw_electrode_names[i]] = i
end

export Potential

struct Potential
    electrodes::Int
    nx::Int
    ny::Int
    nz::Int
    stride::NTuple{3,Float64}
    origin::NTuple{3,Float64}
    data::Array{Float64,4}
    electrode_index::Dict{String,Int}
    electrode_names::Vector{Vector{String}}
end

function Potential(raw::RawPotential,
                   electrode_names::AbstractVector)
    @assert raw.electrodes == length(raw_electrode_names)
    new_electrodes = length(electrode_names)
    data = Array{Float64}(undef, raw.nz, raw.ny, raw.nx, new_electrodes)
    electrode_index = Dict{String,Int}()
    for i in 1:new_electrodes
        electrodes = electrode_names[i]
        first = true
        for elec in electrodes
            electrode_index[elec] = i
            raw_idx = raw_electrode_index[elec]
            if first
                data[:, :, :, i] .= @view raw.data[:, :, :, raw_idx]
                first = false
            else
                data[:, :, :, i] .+= @view raw.data[:, :, :, raw_idx]
            end
        end
        @assert !first
    end
    return Potential(new_electrodes, raw.nx, raw.ny, raw.nz,
                     raw.stride, raw.origin,
                     data, electrode_index, electrode_names)
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        export $(Symbol(name, "_index_to_axis"))
        $(Symbol(name, "_index_to_axis"))(sol::Potential, i) = (i - 1) * sol.stride[$i] + sol.origin[$i]
        export $(Symbol(name, "_axis_to_index"))
        $(Symbol(name, "_axis_to_index"))(sol::Potential, a) = (a - sol.origin[$i]) / sol.stride[$i] + 1
    end
end

function _aliases_to_names(aliases::AbstractDict{Int,Int})
    # Compute the mapping between id's
    raw_electrodes = length(raw_electrode_names)
    id_map = zeros(Int, raw_electrodes)
    id = 0
    new_electrodes = raw_electrodes - length(aliases)
    electrode_names = Vector{Vector{String}}(undef, new_electrodes)
    for i in 1:raw_electrodes
        if i in keys(aliases)
            continue
        end
        id += 1
        id_map[i] = id
        electrode_names[id] = [raw_electrode_names[i]]
    end
    @assert new_electrodes == id
    for (k, v) in aliases
        # The user should connect directly to the final one
        @assert !(v in keys(aliases))
        id = id_map[v]
        @assert id != 0
        name = raw_electrode_names[k]
        push!(electrode_names[id], name)
    end
    return electrode_names
end

function _aliases_to_names(aliases::AbstractDict{S1,S2} where {S1<:AbstractString,S2<:AbstractString})
    return _aliases_to_names(Dict(raw_electrode_index[k]=>raw_electrode_index[v]
                                  for (k, v) in aliases))
end

function _get_electrode_names(aliases, electrode_names)
    if electrode_names !== nothing
        @assert aliases === nothing
        return electrode_names
    end
    if aliases === nothing
        return [[name] for name in raw_electrode_names]
    end
    return _aliases_to_names(aliases)
end

function import_pillbox_v0(filename; aliases=nothing, electrode_names=nothing)
    return Potential(import_pillbox_v0_raw(filename),
                     _get_electrode_names(aliases, electrode_names))
end

function import_pillbox_v1(filename; aliases=nothing, electrode_names=nothing)
    return Potential(import_pillbox_v1_raw(filename),
                     _get_electrode_names(aliases, electrode_names))
end

function import_pillbox_64(filename; aliases=nothing, electrode_names=nothing)
    return Potential(import_pillbox_64_raw(filename),
                     _get_electrode_names(aliases, electrode_names))
end

const _subarray_T = typeof(@view zeros(0, 0, 0, 1)[:, :, :, 1])

struct FitCache
    fitter::Fitting.PolyFitter{3}
    solution::Potential
    cache::Vector{Fitting.PolyFitCache{3,_subarray_T}}
    function FitCache(fitter::Fitting.PolyFitter{3}, solution::Potential)
        return new(fitter, solution,
                   Vector{Fitting.PolyFitCache{3,_subarray_T}}(undef, solution.electrodes))
    end
end

function Base.get(cache::FitCache, idx::Integer)
    if isassigned(cache.cache, idx)
        return cache.cache[idx]
    end
    fit_cache = Fitting.PolyFitCache(cache.fitter, @view cache.solution.data[:, :, :, idx])
    cache.cache[idx] = fit_cache
    return fit_cache
end

Base.get(cache::FitCache, name::AbstractString) =
    get(cache, cache.solution.electrode_index[name])

Base.get(cache::FitCache, electrode::Union{AbstractString,Integer}, pos::NTuple{3};
         fit_center=pos) =
             get(get(cache, electrode), pos; fit_center=fit_center)

get_single(cache::FitCache, electrode::Union{AbstractString,Integer}, pos::NTuple{3},
           orders::NTuple{3}; fit_center=pos) =
               get_single(get(cache, electrode), pos, orders; fit_center=fit_center)

function get_multi_electrodes(cache::FitCache, electrodes_voltages, pos::NTuple{3})
    local res
    for (ele, v) in electrodes_voltages
        term = get(cache, ele, pos) * v
        if !@isdefined(res)
            res = term
        else
            res += term
        end
    end
    if !@isdefined(res)
        return Fitting.PolyFitResult{3}(cache.fitter.orders,
                                        zeros(prod(cache.fitter.orders .+ 1)))
    end
    return res
end

function get_multi_electrodes(cache::FitCache, electrodes_voltages, pos::NTuple{3},
                              orders::NTuple{3})
    res = 0.0
    for (ele, v) in electrodes_voltages
        res += get_single(cache, ele, pos, orders) * v
    end
    return res
end

end
