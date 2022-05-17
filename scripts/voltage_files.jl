#!/usr/bin/julia

function with_write(cb, file::IO)
    cb(file)
end

function with_write(cb, file::AbstractString)
    open(cb, file, "w")
end

function str_float(v)
    if isinteger(v)
        return string(Int(v))
    end
    return uppercase(string(v))
end

struct MapFile
    # Only the names seem to be useful
    names::Vector{String}
    MapFile(names=String[]) = new(names)
end

function load_file(file, ::Type{MapFile})
    res = MapFile()
    for line in eachline(file)
        push!(res.names, split(line, '\t', limit=2)[1])
    end
    return res
end

function write_file(file, data::MapFile)
    with_write(file) do fh
        id = 0
        for name in data.names
            print(fh, "$(name)\t$(id)\t$(id)\r\n")
            id += 1
        end
    end
end

struct CompensationFile
    map::MapFile
    term_names::Vector{String}
    term_values::Vector{Vector{Float64}}
end

function load_file(file, ::Type{CompensationFile}; mapfile=nothing)
    # Assume the eachline iterator is stateful
    term_names = String[]
    lines = eachline(file)
    for line in lines
        fields = split(line, '\t')
        assign_expr = split(fields[1], '=', limit=2)
        if length(assign_expr) < 2
            # We found the last line of the term names
            # and now it's the line with the electrode names
            if mapfile === nothing
                mapfile = MapFile(fields)
            else
                @assert mapfile.names == fields
            end
            break
        end
        # Read term names
        push!(term_names, assign_expr[1])
    end
    mapfile = mapfile::MapFile
    # Read the rest of the lines
    term_values = Vector{Float64}[]
    for line in lines
        fields = parse.(Float64, split(line, '\t'))
        @assert length(fields) == length(mapfile.names)
        push!(term_values, fields)
    end
    @assert length(term_values) == length(term_names)
    return CompensationFile(mapfile, term_names, term_values)
end

function write_file(file, data::CompensationFile)
    with_write(file) do fh
        names_suffix = '\t' ^ (length(data.map.names) - 1)
        id = 0
        for name in data.term_names
            # \n line end seems to work fine with the compensation file
            print(fh, "$(name)=$(id)$(names_suffix)\n")
            id += 1
        end
        println(fh, join(data.map.names, '\t'))
        for values in data.term_values
            println(fh, join((str_float(val) for val in values), '\t'))
        end
    end
end
