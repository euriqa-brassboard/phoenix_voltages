#!/usr/bin/julia

module Potentials

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

function Potential(raw::RawPotential, aliases::Dict{Int,Int})
    # Compute the mapping between id's
    id_map = zeros(Int, raw.electrodes)
    id = 0
    new_electrodes = raw.electrodes - length(aliases)
    data = Array{Float64}(undef, raw.nz, raw.ny, raw.nx, new_electrodes)
    electrode_names = Vector{Vector{String}}(undef, new_electrodes)
    electrode_index = Dict{String,Int}()
    for i in 1:raw.electrodes
        if i in keys(aliases)
            continue
        end
        id += 1
        id_map[i] = id
        data[:, :, :, id] .= @view raw.data[:, :, :, i]
        name = raw_electrode_names[i]
        electrode_index[name] = id
        electrode_names[id] = [name]
    end
    @assert new_electrodes == id
    for (k, v) in aliases
        # The user should connect directly to the final one
        @assert !(v in keys(aliases))
        id = id_map[v]
        @assert id != 0
        data[:, :, :, id] .= @view(data[:, :, :, id]) .+ @view(raw.data[:, :, :, k])
        name = raw_electrode_names[k]
        electrode_index[name] = id
        push!(electrode_names[id], name)
    end
    return Potential(new_electrodes, raw.nx, raw.ny, raw.nz,
                     raw.stride, raw.origin,
                     data, electrode_index, electrode_names)
end

function Potential(raw::RawPotential, aliases::Dict{String,String})
    return Potential(raw, Dict(raw_electrode_index[k]=>raw_electrode_index[v]
                               for (k, v) in aliases))
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        export $(Symbol(name, "_index_to_axis"))
        $(Symbol(name, "_index_to_axis"))(sol::Potential, i) = (i - 1) * sol.stride[$i] + sol.origin[$i]
        export $(Symbol(name, "_axis_to_index"))
        $(Symbol(name, "_axis_to_index"))(sol::Potential, a) = (a - sol.origin[$i]) / sol.stride[$i] + 1
    end
end

function import_pillbox_v0(filename; aliases=Dict{Int,Int}())
    return Potential(import_pillbox_v0_raw(filename), aliases)
end

function import_pillbox_v1(filename; aliases=Dict{Int,Int}())
    return Potential(import_pillbox_v1_raw(filename), aliases)
end

function import_pillbox_64(filename; aliases=Dict{Int,Int}())
    return Potential(import_pillbox_64_raw(filename), aliases)
end

end
