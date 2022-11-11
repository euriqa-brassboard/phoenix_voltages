#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using PyPlot
using MAT

const coeff_data = matopen(joinpath(@__DIR__, "../data/merge_coeff_20221021.mat")) do mat
    return read(mat, "data")
end
const nsection = length(coeff_data)

const transfer = matopen(joinpath(@__DIR__, "../data/merge_20221021_blk70.mat")) do mat
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

const xpos_ums1 = -3045:105
const loading_pos_um = -3045
@assert length(transfer.transfer_solutions[1]) == length(xpos_ums1)

function pos_to_name(xpos_um)
    if xpos_um == loading_pos_um
        return "Loading"
    elseif xpos_um == 0
        return "Center" # Use the hardcoded center position for now.
    else
        return "$(xpos_um)um"
    end
end

function generate_lines_group1(out_lines, xml_io, lines)
    local prev_pos
    local prev_lineidx
    found_center = false
    for (xpos_um, line) in zip(xpos_ums1, lines)
        push!(out_lines, line)
        if xpos_um % 15 != 0
            continue
        end
        lineidx = length(out_lines) - 1
        if xpos_um == 0
            found_center = true
        end
        if @isdefined(prev_pos)
            name = pos_to_name(prev_pos)
            name2 = pos_to_name(xpos_um)
            println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(prev_lineidx)\" startName=\"$(name)\" stopLine=\"$(lineidx)\" stopName=\"$(name2)\"/>")
        end
        prev_pos = xpos_um
        prev_lineidx = lineidx
    end
    @assert found_center
    return
end

function generate_lines(xml_io, transfer_lines)
    out_lines = Vector{Float64}[]
    # XML header
    print(xml_io, """
<?xml version="1.0" ?>
<VoltageAdjust>
  <ShuttlingGraph>
""")
    generate_lines_group1(out_lines, xml_io, transfer_lines[1])
    loading_line = out_lines[1]
    section1_end = length(out_lines) - 1

    load_end_name = pos_to_name(xpos_ums1[end])
    merge_end_name = "200um"
    move_end_name = "300um"

    for line in transfer_lines[2]
        push!(out_lines, line)
    end
    push!(out_lines, transfer_lines[3][end])
    section2_end = length(out_lines) - 1
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(section1_end)\" startName=\"$(load_end_name)\" stopLine=\"$(section2_end)\" stopName=\"$(merge_end_name)\"/>")

    section3_start = length(out_lines)
    push!(out_lines, transfer_lines[1][end])
    for line in transfer_lines[3]
        push!(out_lines, line)
    end
    section3_end = length(out_lines) - 1
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(section3_start)\" startName=\"$(load_end_name)\" stopLine=\"$(section3_end)\" stopName=\"$(move_end_name)\"/>")

    for line in out_lines
        for v in line
            if v > 10 || v < -10
                @warn "Voltage out of range: $v"
            end
        end
    end

    npos = length(out_lines)
    push!(out_lines, loading_line) # LoadingJump
    push!(out_lines, loading_line) # Loading
    push!(out_lines, zeros(length(mapfile.names))) # Zeros
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos)\" startName=\"LoadingJump\" stopLine=\"$(npos + 1)\" stopName=\"Loading\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 1)\" startName=\"Loading\" stopLine=\"$(npos + 2)\" stopName=\"Zeros\"/>")
    print(xml_io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
    return out_lines
end

function generate_output(prefix, scale)
    transfer_lines = [[get_transfer_line(data, scale) for data in solutions]
                      for solutions in transfer.transfer_solutions]
    transfer_file = TransferFile(mapfile, open(io->generate_lines(io, transfer_lines),
                                               "$(prefix).xml", "w"))
    write_file("$(prefix).txt", transfer_file)
    return
end

generate_output("$(outputprefix)_0.2", 0.2)
generate_output("$(outputprefix)_0.1", 0.1)
generate_output("$(outputprefix)_0.05", 0.05)
generate_output("$(outputprefix)_0.025", 0.025)
