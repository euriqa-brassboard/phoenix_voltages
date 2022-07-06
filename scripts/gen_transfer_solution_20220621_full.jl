#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using NaCsPlot
using PyPlot
using MAT

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const short_map = Solutions.load_short_map(
    joinpath(@__DIR__, "../data/electrode_short_202206.csv"))

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file, aliases=short_map)
const fits_cache = Solutions.compensate_fitter1_2(solution)

const mapfile = load_file(ARGS[2], MapFile)
const outputprefix = ARGS[3]

function get_rf_center(xpos_um)
    xidx = Solutions.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function solve_compensate(xpos_um)
    @show xpos_um
    return Solutions.solve_compensate1_nozx(fits_cache, get_rf_center(xpos_um))
end

function get_transfer_line(term, s)
    data = (term[2].x2 .+ term[2].yz) .* s
    return Solutions.get_data_line(solution, mapfile, term[1], data)
end

const xpos_ums = -3080:1155
const loading_pos_um = -3045
const scales = fill(0.25, length(xpos_ums))
scales[1:36] = range(0.2, 0.25, 36)
const comp_terms = [solve_compensate(xpos_um) for xpos_um in xpos_ums]
const lines = [get_transfer_line(term, s) for (term, s) in zip(comp_terms, scales)]

function pos_to_name(xpos_um)
    if xpos_um == loading_pos_um
        return "Loading"
    elseif xpos_um == 0
        return "Center"
    else
        return "$(xpos_um)um"
    end
end

function generate_xml!(io)
    local loading_line
    found_center = false
    # XML header
    print(io, """
<?xml version="1.0" ?>
<VoltageAdjust>
  <ShuttlingGraph>
""")
    npos = length(xpos_ums)
    for i in 1:npos - 1
        xpos_um = Int(xpos_ums[i])
        name = pos_to_name(xpos_um)
        name2 = pos_to_name(xpos_ums[i + 1])
        if xpos_um == loading_pos_um
            loading_line = lines[i]
        elseif xpos_um == 0
            found_center = true
        end
        println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"10\" startLength=\"0\" stopLength=\"0\" startLine=\"$(i - 1)\" startName=\"$(name)\" stopLine=\"$(i)\" stopName=\"$(name2)\"/>")
    end
    @assert @isdefined(loading_line) && found_center
    push!(lines, loading_line) # LoadingJump
    push!(lines, loading_line) # Loading
    push!(lines, zeros(length(mapfile.names))) # Zeros
    println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos)\" startName=\"LoadingJump\" stopLine=\"$(npos + 1)\" stopName=\"Loading\"/>")
    println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"$(npos + 1)\" startName=\"Loading\" stopLine=\"$(npos + 2)\" stopName=\"Zeros\"/>")
    print(io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
end

open(generate_xml!, "$(outputprefix).xml", "w")

const transfer_file = TransferFile(mapfile, lines)

write_file("$(outputprefix).txt", transfer_file)
