#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT
import Ipopt
using JuMP

const coeff_data_nozx, electrode_names_nozx = matopen(joinpath(@__DIR__, "../data/compensate_nozx_20220920.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const coeff_data_zx, electrode_names_zx = matopen(joinpath(@__DIR__, "../data/compensate_20220920.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

@assert size(coeff_data_zx) == size(coeff_data_nozx)
@assert electrode_names_zx == electrode_names_nozx
@assert [data["xpos_um"] for data in coeff_data_zx] == [data["xpos_um"] for data in coeff_data_nozx]

# Stitch the zx vs no-zx solutions together at around -950 um,
# where the effect of relaxing the zx term should be minimal
# at least on the X2 anx YZ terms.
# The global smoothing should hopefully make any transition in the crossover area
# smooth.
const switchover_idx = findfirst(x->x["xpos_um"] >= -950, coeff_data_zx)

const coeff_data = [coeff_data_nozx[1:switchover_idx - 1];
                    coeff_data_zx[switchover_idx:end]]
const electrode_names = electrode_names_nozx

@assert size(coeff_data) == size(coeff_data_nozx)
@assert [data["xpos_um"] for data in coeff_data] == [data["xpos_um"] for data in coeff_data_nozx]

const prefix = joinpath(@__DIR__, "../data/transfer_smooth_20220920")

function solve_all(diff_weight, diff_flatten_weight)
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
        maxmaxweight = 0
    end
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
            if diff_flatten_weight > 0 && -3050 < data["xpos_um"] < 1130
                @constraint(model, maxmaxdiff >= maxdiff)
                maxmaxweight += 1
            end
        end
        prev_vmap = vmap
    end
    if diff_flatten_weight > 0
        @objective(model, Min, sum(maxvs) + sum(maxdiffs) * diff_weight + maxmaxdiff * (maxmaxweight * diff_flatten_weight))
    else
        @objective(model, Min, sum(maxvs) + sum(maxdiffs) * diff_weight)
    end
    JuMP.optimize!(model)
    return [[value(v) for v in x] for x in xs]
end

function pack_data(data, vals)
    solution = data["solution"]
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end

matopen("$(prefix).mat", "w") do mat
    voltages = @time(solve_all(0.05, 1))
    transfer_solutions = [pack_data(data, vals) for (data, vals)
                              in zip(coeff_data, voltages)]
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end

matopen("$(prefix)_noglobal.mat", "w") do mat
    voltages = @time(solve_all(0.05, 0))
    transfer_solutions = [pack_data(data, vals) for (data, vals)
                              in zip(coeff_data, voltages)]
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
