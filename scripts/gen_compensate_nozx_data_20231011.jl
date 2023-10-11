#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT
import Ipopt
using JuMP

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/compensate_nozx_20220927.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const outputdir = joinpath(@__DIR__, "../data/compensate_nozx_20220927")

function solve_all(termidx)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", 30000)
    set_optimizer_attribute(model, "print_level", 5)
    xs = Vector{AffExpr}[]
    maxvs = VariableRef[]
    for data in coeff_data
        solution = data["solution"]
        x0 = solution[:, termidx]
        B = data["free_solution"]
        nx, nt = size(B)
        @assert nx == length(x0)
        t = @variable(model, [1:nt])
        x = @expression(model, B * t .+ x0)
        maxv = @variable(model)
        @constraint(model, maxv .>= x)
        @constraint(model, maxv .>= .-x)
        push!(xs, x)
        push!(maxvs, maxv)
    end
    @objective(model, Min, sum(maxvs))
    JuMP.optimize!(model)
    return [[value(v) for v in x] for x in xs]
end

const term_names = ("dx", "dy", "dz", "xy", "yz", "z2", "x2", "x3", "x4")

function pack_data(data, vals)
    solution = data["solution"]
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end

function gen_solution(termidx)
    mkpath(outputdir)
    matopen(joinpath(outputdir, "$(termidx).mat"), "w") do mat
        voltages = @time(solve_all(termidx))
        transfer_solutions = [pack_data(data, vals) for (data, vals)
                                  in zip(coeff_data, voltages)]
        write(mat, "electrode_names", electrode_names)
        write(mat, "transfer_solutions", transfer_solutions)
        write(mat, "termidx", termidx)
        write(mat, "termname", term_names[termidx])
    end
end

gen_solution(parse(Int, ARGS[1]))
