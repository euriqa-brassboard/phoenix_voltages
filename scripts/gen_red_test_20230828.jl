#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
using PhoenixVoltages.Outputs

const data_prefix = joinpath(@__DIR__, "../data/red_chamber_20230828/")
const outputprefix = joinpath(data_prefix, "test_brassboard")

const red_trans = load_file(joinpath(data_prefix, "Phoenix_2_m35.txt"), TransferFile)
const red_comp = load_file(joinpath(data_prefix, "global_adjust_2_m35.txt"),
                           CompensationFile)

const mapfile = load_file(ARGS[1], MapFile)

function _map_line(line, old_index, names)
    return [begin
                idx = get(old_index, name, 0)
                idx == 0 ? 0.0 : line[idx]
            end for name in names]
end

function map_line(line, old_map, new_map)
    old_index = Dict(name=>i for (i, name) in enumerate(old_map.names))
    return _map_line(line, old_index, new_map.names)
end

function map_lines(lines, old_map, new_map)
    old_index = Dict(name=>i for (i, name) in enumerate(old_map.names))
    return [_map_line(line, old_index, new_map.names) for line in lines]
end

const test_comp_terms = Dict(name=>map_line(line, red_comp.map, mapfile)
                            for (name, line)
                                in zip(red_comp.term_names, red_comp.term_values))
const test_trans_lines = map_lines(red_trans.line_values, red_trans.map, mapfile)

const test_line1 = test_trans_lines[1] .+ 0.1 .* test_comp_terms["simpleRot"]
const test_line2 = test_trans_lines[1] .+ 0.3 .* test_comp_terms["simpleRot"]
const test_line3 = test_trans_lines[1] .+ 0.5 .* test_comp_terms["simpleRot"]
const test_line4 = test_trans_lines[1] .+ 1.0 .* test_comp_terms["simpleRot"]
const test_line0 = test_trans_lines[1]

function generate_lines(xml_io)
    # XML header
    lines = Vector{Float64}[]
    print(xml_io, """
<?xml version="1.0" ?>
<VoltageAdjust>
  <ShuttlingGraph>
""")
    push!(lines, test_line1)
    push!(lines, test_line2)
    push!(lines, test_line3)
    push!(lines, test_line4)
    push!(lines, test_line0)
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"0\" startName=\"LoadingJump\" stopLine=\"1\" stopName=\"Loading\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"1\" startName=\"Loading\" stopLine=\"2\" stopName=\"Center\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"2\" startName=\"Center\" stopLine=\"3\" stopName=\"Center2\"/>")
    println(xml_io, "    <ShuttleEdge idleCount=\"500\" steps=\"3\" startLength=\"0\" stopLength=\"0\" startLine=\"3\" startName=\"Center2\" stopLine=\"4\" stopName=\"Zeros\"/>")
    print(xml_io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
    return lines
end

const transfer_file = TransferFile(mapfile, open(generate_lines,
                                                 "$(outputprefix).xml", "w"))

write_file("$(outputprefix).txt", transfer_file)
