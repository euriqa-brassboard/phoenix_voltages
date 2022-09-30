#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT
import Ipopt
using JuMP

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/merge_coeff_20220929.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end
const xpos_ums = [data["xpos_um"] for data in coeff_data]

const prefix = joinpath(@__DIR__, "../data/merge_smooth_20220929")

function solve_all(diff_weight, diff_flatten_weight, max_flatten_weight)
    @show diff_weight
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", 30000)
    set_optimizer_attribute(model, "print_level", 5)
    xs = Vector{AffExpr}[]
    maxvs = VariableRef[]
    maxdiffs = VariableRef[]
    local prev_vmap
    if diff_flatten_weight > 0
        @variable(model, maxmaxdiff)
        maxmaxdiffweight = 0
    end
    if max_flatten_weight > 0
        @variable(model, maxmaxv)
        maxmaxvweight = 0
    end
    for data in coeff_data
        x0 = data["solution"]
        B = data["free_solution"]
        nx, nt = size(B)
        @assert nx == length(x0)
        t = @variable(model, [1:nt])
        if "limited_solution" in keys(data)
            limited = data["limited_solution"]
            limited = reshape(limited, size(limited, 1), :)
            @assert nx == size(limited, 1)
            nl = size(limited, 2)
            l = @variable(model, [1:nl])
            @constraint(model, l .<= 1)
            @constraint(model, l .>= -1)
            x = @expression(model, B * t .+ x0 .+ limited * l)
        else
            x = @expression(model, B * t .+ x0)
        end
        maxv = @variable(model)
        @constraint(model, maxv .>= x)
        @constraint(model, maxv .>= .-x)
        push!(xs, x)
        push!(maxvs, maxv)
        vmap = Dict(zip(data["electrodes"], x))
        if max_flatten_weight > 0 && -2500 < data["xpos_um"] < 1000
            @constraint(model, maxmaxv >= maxv)
            maxmaxvweight += 1
        end
        if @isdefined(prev_vmap)
            maxdiff = @variable(model)
            for k in union(keys(vmap), keys(prev_vmap))
                v1 = get(vmap, k, 0)
                v2 = get(prev_vmap, k, 0)
                @constraint(model, maxdiff >= v1 - v2)
                @constraint(model, maxdiff >= v2 - v1)
            end
            push!(maxdiffs, maxdiff)
            if diff_flatten_weight > 0
                @constraint(model, maxmaxdiff >= maxdiff)
                maxmaxdiffweight += 1
            end
        end
        prev_vmap = vmap
    end
    obj = sum(maxvs) + sum(maxdiffs) * diff_weight
    if diff_flatten_weight > 0
        obj = obj + maxmaxdiff * (maxmaxdiffweight * diff_flatten_weight)
    end
    if max_flatten_weight > 0
        obj = obj + maxmaxv * (maxmaxvweight * max_flatten_weight)
    end
    @objective(model, Min, obj)
    JuMP.optimize!(model)
    return [[value(v) for v in x] for x in xs]
end

function pack_data(data, vals)
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end

function gen_global()
    voltages = @time(solve_all(0.05, 2, 10))
    transfer_solutions = [pack_data(data, vals) for (data, vals)
                              in zip(coeff_data, voltages)]
    matopen("$(prefix).mat", "w") do mat
        write(mat, "electrode_names", electrode_names)
        write(mat, "transfer_solutions", transfer_solutions)
    end
end

function gen_noglobal()
    voltages = @time(solve_all(0.05, 0, 10))
    transfer_solutions = [pack_data(data, vals) for (data, vals)
                              in zip(coeff_data, voltages)]
    matopen("$(prefix)_noglobal.mat", "w") do mat
        write(mat, "electrode_names", electrode_names)
        write(mat, "transfer_solutions", transfer_solutions)
    end
end

if isempty(ARGS)
    gen_global()
    gen_noglobal()
elseif ARGS[1] == "g"
    gen_global()
elseif ARGS[1] == "ng"
    gen_noglobal()
else
    error("Invalid argument: $(ARGS)")
end
