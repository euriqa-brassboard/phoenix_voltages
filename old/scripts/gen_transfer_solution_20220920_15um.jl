#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end

const transfer = matopen(joinpath(@__DIR__, "../data/transfer_smooth_20220920.mat")) do mat
    return (transfer_solutions=read(mat, "transfer_solutions"),
            electrode_names=read(mat, "electrode_names"))
end
const electrode_names = Vector{String}.(transfer.electrode_names)

const mapfile = load_file(ARGS[1], MapFile)
const outputprefix = ARGS[2]

function get_transfer_line(data, scale=1.0)
    voltage_map = Dict{String,Float64}()
    for (idx, v) in zip(data["electrodes"], data["voltages"])
        # This is the factor that makes sure everything is within +-10
        # for the current solution.
        v = v * scale
        for name in electrode_names[idx]
            voltage_map[name] = v
        end
    end
    nelectrodes = length(mapfile.names)
    values = Vector{Float64}(undef, nelectrodes)
    for i in 1:nelectrodes
        values[i] = get(voltage_map, mapfile.names[i], 0.0)
    end
    return values
end

# 0.5 here corresponds to a X2 of 0.5 and QYZ of 0.375.
const transfer_lines = [(data["xpos_um"], get_transfer_line(data, 0.5))
                        for data in transfer.transfer_solutions]

const xpos_ums = -3075:15:500
const loading_pos_um = -3045

function pos_to_name(xpos_um)
    if xpos_um == loading_pos_um
        return "Loading"
    elseif xpos_um == 0
        return "Center"
    else
        return "$(xpos_um)um"
    end
end

function generate_lines(xml_io, transfer_lines)
    node_idx = 1
    node_xpos = xpos_ums[node_idx]
    lines = Vector{Float64}[]
    local loading_line
    local prev_lineidx
    found_center = false
    # XML header
    print(xml_io, """
<?xml version="1.0" ?>
<VoltageAdjust>
  <ShuttlingGraph>
""")
    for line in transfer_lines
        if line[1] < node_xpos
            if node_idx > 1
                push!(lines, line[2])
            end
            continue
        end
        push!(lines, line[2])
        @assert line[1] == node_xpos
        lineidx = length(lines) - 1
        if node_idx > 1
            name = pos_to_name(xpos_ums[node_idx - 1])
            name2 = pos_to_name(node_xpos)
            if node_xpos == loading_pos_um
                loading_line = line[2]
            elseif node_xpos == 0
                found_center = true
            end
            println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(prev_lineidx)\" startName=\"$(name)\" stopLine=\"$(lineidx)\" stopName=\"$(name2)\"/>")
        end
        prev_lineidx = lineidx
        node_idx += 1
        if node_idx > length(xpos_ums)
            break
        end
        node_xpos = xpos_ums[node_idx]
    end

    for line in lines
        for v in line
            if v > 10 || v < -10
                @warn "Voltage out of range: $v"
            end
        end
    end

    @assert @isdefined(loading_line) && found_center
    npos = length(lines)
    push!(lines, loading_line) # LoadingJump
    push!(lines, loading_line) # Loading
    push!(lines, zeros(length(mapfile.names))) # Zeros
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos)\" startName=\"LoadingJump\" stopLine=\"$(npos + 1)\" stopName=\"Loading\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 1)\" startName=\"Loading\" stopLine=\"$(npos + 2)\" stopName=\"Zeros\"/>")
    print(xml_io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
    return lines
end

const transfer_file = TransferFile(mapfile, open(io->generate_lines(io, transfer_lines),
                                                 "$(outputprefix).xml", "w"))

write_file("$(outputprefix).txt", transfer_file)
