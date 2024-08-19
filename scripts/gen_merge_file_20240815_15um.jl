#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using PyPlot
using MAT

const centers = Solutions.CenterTracker()

const transfer = matopen(joinpath(@__DIR__, "../data/merge_20240815_blk70.mat")) do mat
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

const xpos_ums = -3045:15:0
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

function gen_test_line(names)
    nelectrodes = length(names)
    line = zeros(nelectrodes)
    i1 = 1
    i2 = nelectrodes
    v1 = -99
    v2 = 99
    while i1 < i2
        if !startswith(names[i1], "GND")
            line[i1] = v1 / 10
            v1 += 2
        end
        if !startswith(names[i2], "GND")
            line[i2] = v2 / 10
            v2 -= 2
        end
        i1 += 1
        i2 -= 1
    end
    if i1 == i2
        if !startswith(names[i1], "GND")
            line[i1] = v1 / 10
            v1 += 2
        end
    end
    if v1 >= v2
        @warn "None unique voltages on the test line"
    end
    return line
end

function gen_repel_line(names, v)
    nelectrodes = length(names)
    line = zeros(nelectrodes)
    for i in 1:nelectrodes
        if !startswith(names[i], "GND")
            line[i] = v
        end
    end
    return line
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
        if node_xpos == loading_pos_um
            loading_line = line[2]
        elseif node_xpos == 0
            found_center = true
        end
        if node_idx > 1
            name = pos_to_name(xpos_ums[node_idx - 1])
            name2 = pos_to_name(node_xpos)
            println(xml_io, "    <ShuttleEdge idleCount=\"5000\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"$(prev_lineidx)\" startName=\"$(name)\" stopLine=\"$(lineidx)\" stopName=\"$(name2)\"/>")
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
    push!(lines, gen_test_line(mapfile.names)) # Test
    push!(lines, gen_repel_line(mapfile.names, 9.9)) # RepelP
    push!(lines, gen_repel_line(mapfile.names, -9.9)) # RepelM
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos)\" startName=\"LoadingJump\" stopLine=\"$(npos + 1)\" stopName=\"Loading\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 1)\" startName=\"Loading\" stopLine=\"$(npos + 2)\" stopName=\"Zeros\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 2)\" startName=\"Zeros\" stopLine=\"$(npos + 3)\" stopName=\"DCTest\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 3)\" startName=\"DCTest\" stopLine=\"$(npos + 4)\" stopName=\"Repel\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 4)\" startName=\"Repel\" stopLine=\"$(npos + 5)\" stopName=\"RepelM\"/>")
    print(xml_io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
    return lines
end

function generate_output(prefix, scale)
    transfer_lines = [(data["xpos_um"], get_transfer_line(data, scale))
                      for data in transfer.transfer_solutions]
    transfer_file = TransferFile(mapfile, open(io->generate_lines(io, transfer_lines),
                                               "$(prefix).xml", "w"))
    write_file("$(prefix).txt", transfer_file)
    return
end

# generate_output("$(outputprefix)_0.08", 0.08)
generate_output("$(outputprefix)_0.05", 0.05)
generate_output("$(outputprefix)_0.025", 0.025)
