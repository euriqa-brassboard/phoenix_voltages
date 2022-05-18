#!/usr/bin/julia

struct VoltageSolutions
    electrodes::Int
    nx::Int
    ny::Int
    nz::Int
    vsets::Int
    xaxis::NTuple{3,Float64}
    yaxis::NTuple{3,Float64}
    stride::NTuple{3,Float64}
    origin::NTuple{3,Float64}
    elecrodemapping::Vector{Int}
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
function import_pillbox_v0(filename)
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
        return VoltageSolutions(electrodes, nx, ny, nz, vsets,
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
function import_pillbox_v1(filename)
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
        return VoltageSolutions(electrodes, nx, ny, nz, vsets, xaxis, yaxis,
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
function import_pillbox_64(filename)
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
        return VoltageSolutions(electrodes, nx, ny, nz, vsets, xaxis, yaxis,
                                stride, origin, electrodemapping, data)
    end
end

for (name, i) in ((:x, 1), (:y, 2), (:z, 3))
    @eval begin
        $(Symbol(name, "_index_to_axis"))(sol::VoltageSolutions, i) = (i - 1) * sol.stride[$i] + sol.origin[$i]
        $(Symbol(name, "_axis_to_index"))(sol::VoltageSolutions, a) = (a - sol.origin[$i]) / sol.stride[$i] + 1
    end
end
