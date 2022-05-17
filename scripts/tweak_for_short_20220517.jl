#!/usr/bin/julia

include("voltage_files.jl")

# Electrode index eidx1 and eidx2 are shorted together and must have the same voltage
# For line index <= lidx1, we'll use the voltage on eidx1 on both.
# For line index >= lidx2, we'll use the voltage on eidx2 on both.
# For lidx1 < line index < lidx2, we'll use a linear interpolation
function select_closer_electrode!(data::TransferFile, eidx1, eidx2, lidx1, lidx2)
    nlines = length(data.line_values)
    for lineno in 1:nlines
        line_value = data.line_values[lineno]
        v1 = line_value[eidx1]
        v2 = line_value[eidx2]
        if lineno <= lidx1
            v = v1
        elseif lineno >= lidx2
            v = v2
        else
            v = (v1 * (lidx2 - lineno) + v2 * (lineno - lidx1)) / (lidx2 - lidx1)
        end
        line_value[eidx1] = v
        line_value[eidx2] = v
    end
end

function tweak!(data::TransferFile)
    name_map = Dict(data.map)
    select_closer_electrode!(data, name_map["L2"], name_map["L6"],
                             107, 247)
    select_closer_electrode!(data, name_map["L9"], name_map["Q7"],
                             1000, 2000)
end

const file = load_file(ARGS[1], TransferFile)
tweak!(file)
write_file(ARGS[2], file)
