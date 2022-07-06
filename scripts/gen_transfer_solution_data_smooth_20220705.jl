#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT
import Ipopt
using JuMP

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/compensate_nozx_20220705.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const prefix = joinpath(@__DIR__, "../data/transfer_smooth_20220705")

function solve_all(diff_weight)
    @show diff_weight
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_cpu_time", 3600.0)
    set_optimizer_attribute(model, "print_level", 5)
    xs = Vector{AffExpr}[]
    maxvs = VariableRef[]
    maxdiffs = VariableRef[]
    local prev_vmap
    @variable(model, maxmaxdiff)
    maxmaxweight = 0
    for data in coeff_data
        solution = data["solution"]
        x0 = solution[:, 7] .+ solution[:, 5] .* 0.75
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
        vmap = Dict(zip(data["electrodes"], x))
        if @isdefined(prev_vmap)
            maxdiff = @variable(model)
            for k in union(keys(vmap), keys(prev_vmap))
                v1 = get(vmap, k, 0)
                v2 = get(prev_vmap, k, 0)
                @constraint(model, maxdiff >= v1 - v2)
                @constraint(model, maxdiff >= v2 - v1)
            end
            push!(maxdiffs, maxdiff)
            if -3050 < data["xpos_um"] < 1130
                @constraint(model, maxmaxdiff >= maxdiff)
                maxmaxweight += 1
            end
        end
        prev_vmap = vmap
    end
    @objective(model, Min, sum(maxvs) + sum(maxdiffs) * diff_weight + maxmaxdiff * (maxmaxweight))
    JuMP.optimize!(model)
    return [[value(v) for v in x] for x in xs]
end

function pack_data(data, vals)
    solution = data["solution"]
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end

const voltages = @time(solve_all(0.05))
const transfer_solutions = [pack_data(data, vals) for (data, vals)
                                in zip(coeff_data, voltages)]
matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
