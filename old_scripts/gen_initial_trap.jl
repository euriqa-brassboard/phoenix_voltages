#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using Printf
using DataStructures
using PhoenixVoltages.Mappings

const prefix = joinpath(@__DIR__, "../voltages/initial")

function ao2electrode(ao)
    if ao % 25 == 12
        return "GND$(ao ÷ 25)"
    else
        return replace(ao_to_electrode[ao], "-"=>"_")
    end
end

function ao2dsub(ao)
    dsub_name = ("N", "S", "E", "W")[ao ÷ 25 + 1]
    dsub_pin = ao % 25 + 1
    return "$(dsub_name)$(@sprintf("%02d", dsub_pin))"
end

open("$(prefix)_electrode_mapping.txt", "w") do io
    for i in 1:100
        ao = i - 1
        electrode = ao2electrode(ao)
        dsub = ao2dsub(ao)
        println(io, "$electrode\t$ao\t$dsub")
    end
end

open("$(prefix)_voltage_definition.txt", "w") do io
    for i in 1:100
        ao = i - 1
        electrode = ao2electrode(ao)
        if i != 1
            print(io, "\t")
        end
        print(io, electrode)
    end
    println(io)
    for i in 1:100
        if i != 1
            print(io, "\t")
        end
        print(io, "0")
    end
    println(io)
end

const solutions = OrderedDict("M_X2"=>Dict("O0"=>-1.89836, "O1"=>-1.89836,
                                           "L4"=>-0.983939, "L5"=>-0.983939,
                                           "L0"=>11.1849, "L1"=>11.1849,
                                           "L8"=>11.1849, "L9"=>11.1849),
                              "M_DZ"=>Dict("O0"=>-0.777, "O1"=>-0.777),
                              "M_QYZ"=>Dict("O0"=>-0.32994, "O1"=>0.32994,
                                            "L4"=>0.384836, "L5"=>-0.384836,
                                            "L0"=>1.07679, "L8"=>1.07679,
                                            "L1"=>-1.07679, "L9"=>-1.07679))

for electrode in keys(electrode_to_ao)
    solutions["E_$electrode"] = Dict(electrode=>1.0)
end

open("$(prefix)_global_adjust.txt", "w") do io
    names = collect(keys(solutions))
    for i in 1:length(names)
        name = names[i]
        println(io, "$(name)=$(i - 1)")
    end
    for i in 1:100
        ao = i - 1
        electrode = ao2electrode(ao)
        if i != 1
            print(io, "\t")
        end
        print(io, electrode)
    end
    println(io)
    for name in names
        solution = solutions[name]
        for (e, v) in solution
            @assert(e in keys(electrode_to_ao))
        end
        for i in 1:100
            ao = i - 1
            electrode = ao2electrode(ao)
            if i != 1
                print(io, "\t")
            end
            if electrode in keys(solution)
                print(io, solution[electrode])
            else
                print(io, 0)
            end
        end
        println(io)
    end
end
