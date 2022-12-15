#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Outputs
using MAT

const prefixes = ["compensate_red_20221213"=>("origin", 10),
                  "compensate_red2_20221213"=>("fewer", 10)]

const global_adjust_electrodes = [
    "Q00", "Q01", "Q02", "Q03", "Q04", "Q05", "Q06", "Q07", "Q08", "Q09",
    "Q10", "Q11", "Q12", "Q13", "Q14", "Q15", "Q16", "Q17", "Q18", "Q19",
    "Q20", "Q21", "Q22", "Q23", "Q24", "Q25", "Q26", "Q27", "Q28", "Q29",
    "Q30", "Q31", "Q32", "Q33", "Q34", "Q35", "Q36", "Q37", "Q38", "Q39",
    "Q40", "Q41", "Q42", "Q43", "Q44", "Q45", "Q46", "Q47", "Q48", "Q49",
    "Q50", "Q51", "Q52", "Q53", "Q54", "Q55", "Q56", "Q57", "Q58", "Q59",
    "Q60", "Q61", "Q62", "Q63", "Q64", "Q65"]

const trap_electrodes = [
    "L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "O0", "O1",
    "Q00", "Q01", "Q02", "Q03", "Q04", "Q05", "Q06", "Q07", "Q08", "Q09",
    "Q10", "Q11", "Q12", "Q13", "Q14", "Q15", "Q16", "Q17", "Q18", "Q19",
    "Q20", "Q21", "Q22", "Q23", "Q24", "Q25", "Q26", "Q27", "Q28", "Q29",
    "Q30", "Q31", "Q32", "Q33", "Q34", "Q35", "Q36", "Q37", "Q38", "Q39",
    "Q40", "Q41", "Q42", "Q43", "Q44", "Q45", "Q46", "Q47", "Q48", "Q49",
    "Q50", "Q51", "Q52", "Q53", "Q54", "Q55", "Q56", "Q57", "Q58", "Q59",
    "Q60", "Q61", "Q62", "Q63", "Q64", "Q65",
    "S00", "S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11"]

get_data_line(vmap, electrodes) =
    Float64[get(vmap, e, 0.0) for e in electrodes]

function fix_electrode_name(orig_name)
    m = match(r"^Q([0-9])$", orig_name)
    if m !== nothing
        return "Q0$(m[1])"
    end
    m = match(r"^S([0-9])$", orig_name)
    if m !== nothing
        return "S0$(m[1])"
    end
    return orig_name
end

function unfix_electrode_name(orig_name)
    m = match(r"^([A-Z])0([0-9])$", orig_name)
    if m !== nothing
        return "$(m[1])$(m[2])"
    end
    return orig_name
end

function solution_to_vmap(solution, electrode_index, electrode_names)
    vmap = Dict{String,Float64}()
    for (v, ei) in zip(solution, electrode_index)
        for name in electrode_names[ei]
            vmap[fix_electrode_name(name)] = v
        end
    end
    return vmap
end

mutable struct SolutionData
    x_comp::Vector{Float64}
    y_comp::Vector{Float64}
    z_comp::Vector{Float64}
    simpleRot::Vector{Float64}
    PotentialSqueeze::Vector{Float64}
    trap_adjust::Vector{Float64}
    trap::Vector{Float64}
end

function load_solution(name)
    data = matread(joinpath(@__DIR__, "../data", name))
    termname = data["termname"]
    solution = data["transfer_solutions"]
    xpos_um = [sol["xpos_um"] for sol in solution]
    electrode_names = data["electrode_names"]
    vmaps = [solution_to_vmap(sol["voltages"], sol["electrodes"], electrode_names)
             for sol in solution]
    return (termname=termname, xpos_um=xpos_um, vmaps=vmaps)
end

function vmap_to_solution_data(vmaps, idx)
    return SolutionData(
        get_data_line(vmaps["dx"][idx], global_adjust_electrodes) .* 1000,
        get_data_line(vmaps["dy"][idx], global_adjust_electrodes) .* 1000,
        get_data_line(vmaps["dz"][idx], global_adjust_electrodes) .* 1000,
        get_data_line(vmaps["yz"][idx], global_adjust_electrodes),
        get_data_line(vmaps["z2"][idx], global_adjust_electrodes),
        get_data_line(vmaps["x2"][idx], global_adjust_electrodes),
        get_data_line(vmaps["x2"][idx], trap_electrodes) .* 0.09
    )
end

function load_solutions(prefix, num)
    local xpos_um
    vmaps = Dict{String,Vector{Dict{String,Float64}}}()
    for i in 1:num
        sol = load_solution(joinpath(prefix, "$(i).mat"))
        if @isdefined(xpos_um)
            @assert xpos_um == sol.xpos_um
        else
            xpos_um = sol.xpos_um
        end
        vmaps[sol.termname] = sol.vmaps
    end
    sol_datas = [vmap_to_solution_data(vmaps, i) for i in 1:length(xpos_um)]
    return (xpos_um=xpos_um, sol_datas=sol_datas)
end

const solutions = [(name, load_solutions(prefix, num))
                   for (prefix, (name, num)) in prefixes]

const outputprefix = joinpath(@__DIR__, "../data/compensate_red")

const header = readlines(joinpath(@__DIR__, "../data/PhoenixGlobalAdjust_header.txt"))

line_to_txt(line) = join([v == 0 ? "0" : string(v) for v in line], "\t")

function println_w(io, args...)
    print(io, args..., "\r\n")
end

function write_global_adjust(name, info)
    adj_lines = Tuple{String,String}[]
    push!(adj_lines, ("z_comp", line_to_txt(info.z_comp)))
    push!(adj_lines, ("y_comp", line_to_txt(info.y_comp)))
    push!(adj_lines, ("x_comp", line_to_txt(info.x_comp)))
    push!(adj_lines, ("simpleRot", line_to_txt(info.simpleRot)))
    push!(adj_lines, ("PotentialSqueeze", line_to_txt(info.PotentialSqueeze)))
    push!(adj_lines, ("trap", line_to_txt(info.trap_adjust)))
    for (idx, ele) in enumerate(global_adjust_electrodes)
        ele_line = line_to_txt([j == idx ? 1 : 0
                                for j in 1:length(global_adjust_electrodes)])
        push!(adj_lines, (unfix_electrode_name(ele), ele_line))
    end

    open(name, "w") do io
        for line in header
            println_w(io, line)
        end
        for (i, (name, line)) in enumerate(adj_lines)
            println_w(io, "$(name)=$(i - 1)")
        end
        println_w(io)
        println_w(io, join(global_adjust_electrodes, "\t"))
        for (name, line) in adj_lines
            println_w(io, line)
        end
    end
end

function write_trap(name, info)
    open(name, "w") do io
        println_w(io, join(trap_electrodes, "\t"))
        println_w(io, line_to_txt(info.trap))
    end
end

for (name, solution) in solutions
    for i in 1:length(solution.xpos_um)
        outputdir = "$(outputprefix)_$(name)/$(solution.xpos_um[i])um"
        mkpath(outputdir)
        write_global_adjust(joinpath(outputdir, "global_adjust.txt"),
                            solution.sol_datas[i])
        write_trap(joinpath(outputdir, "trap.txt"), solution.sol_datas[i])
    end
end
