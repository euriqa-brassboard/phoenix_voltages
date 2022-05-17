#!/usr/bin/julia

function with_write(cb, file::IO)
    cb(file)
end

function with_write(cb, file::AbstractString)
    open(cb, file, "w")
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
