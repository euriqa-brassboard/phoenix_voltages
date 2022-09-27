#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using MAT

const prefixes = ["compensate_20220927"=>("origin", 10),
                  "compensate_x2z_20220927"=>("x2z", 11),
                  "compensate_nozx_20220927"=>("nozx", 9)]
const mapfile = load_file(ARGS[1], MapFile)

function get_data_line(solution, electrode_index, id_map, mapfile::MapFile)
    nelectrodes = length(mapfile.names)
    values = Vector{Float64}(undef, nelectrodes)
    for i in 1:nelectrodes
        id = get(electrode_index, mapfile.names[i], 0)
        idx = get(id_map, id, 0)
        values[i] = get(solution, idx, 0.0)
    end
    return values
end

function load_solution(name, mapfile::MapFile)
    data = matread(joinpath(@__DIR__, "../data", name))
    termname = data["termname"]
    solution = data["transfer_solutions"]
    electrode_index = Dict{String,Int}()
    for (i, enames) in enumerate(data["electrode_names"])
        for ename in enames
            electrode_index[ename] = i
        end
    end
    xpos_um = [sol["xpos_um"] for sol in solution]
    lines = [get_data_line(sol["voltages"], electrode_index,
                           Dict(id=>idx for (idx, id) in enumerate(sol["electrodes"])),
                           mapfile)
             for sol in solution]
    return (termname=termname, xpos_um=xpos_um, lines=lines)
end

function lines_to_compensation_file(lines, idx, mapfile::MapFile)
    term_names = String[]
    term_values = Vector{Float64}[]
    nelectrodes = length(mapfile.names)
    # For DX, DY, DZ the EURIQA frontend expects a different unit in the config file
    # compared to the UI...
    # DX
    push!(term_names, "DX")
    push!(term_values, lines["dx"][idx] .* 1000)
    # DY
    push!(term_names, "DY")
    push!(term_values, lines["dy"][idx] .* 1000)
    # DZ
    push!(term_names, "DZ")
    push!(term_values, lines["dz"][idx] .* 1000)
    # QZY
    push!(term_names, "QZY")
    push!(term_values, lines["yz"][idx])
    # QZZ
    push!(term_names, "QZZ")
    push!(term_values, lines["z2"][idx])
    # QXZ
    push!(term_names, "QXZ")
    if !haskey(lines, "zx")
        push!(term_values, zeros(nelectrodes))
    else
        push!(term_values, lines["zx"][idx])
    end
    # X1
    # DX is in V/m, X1 is in 525 uV / 2.74 um
    push!(term_names, "X1")
    push!(term_values, lines["dx"][idx] .* (Solutions.V_unit_uV / Solutions.l_unit_um))
    # X2
    push!(term_names, "X2")
    push!(term_values, lines["x2"][idx])
    # X3
    push!(term_names, "X3")
    push!(term_values, lines["x3"][idx])
    # X4
    push!(term_names, "X4")
    push!(term_values, lines["x4"][idx])
    # JunctionCenter
    push!(term_names, "JunctionCenter")
    push!(term_values, zeros(nelectrodes))
    # JunctionTransition
    push!(term_names, "JunctionTransition")
    push!(term_values, zeros(nelectrodes))
    # LoadInners
    push!(term_names, "LoadInners")
    push!(term_values, zeros(nelectrodes))
    # LoadOuters
    push!(term_names, "LoadOuters")
    push!(term_values, zeros(nelectrodes))
    # Eject
    push!(term_names, "Eject")
    push!(term_values, zeros(nelectrodes))

    return CompensationFile(mapfile, term_names, term_values)
end

function load_solutions(prefix, num, mapfile::MapFile)
    local xpos_um
    lines = Dict{String,Vector{Vector{Float64}}}()
    for i in 1:num
        sol = load_solution(joinpath(prefix, "$(i).mat"), mapfile)
        if @isdefined(xpos_um)
            @assert xpos_um == sol.xpos_um
        else
            xpos_um = sol.xpos_um
        end
        lines[sol.termname] = sol.lines
    end
    comp_files = [lines_to_compensation_file(lines, i, mapfile)
                  for i in 1:length(xpos_um)]
    return (xpos_um=xpos_um, comp_files=comp_files)
end

const solutions = [(name, load_solutions(prefix, num, mapfile))
                   for (prefix, (name, num)) in prefixes]

const outputprefix = ARGS[2]

for (name, solution) in solutions
    outputdir = "$(outputprefix)_$(name)"
    mkpath(outputdir)
    for i in 1:length(solution.xpos_um)
        write_file(joinpath(outputdir, "comp_$(solution.xpos_um[i])um.txt"),
                   solution.comp_files[i])
    end
end
