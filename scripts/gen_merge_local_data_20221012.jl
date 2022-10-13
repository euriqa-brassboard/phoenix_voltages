#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using JuMP
using Ipopt
using MAT

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/merge_coeff_20221012.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const xpos_ums = [data["xpos_um"] for data in coeff_data]

const prefix = joinpath(@__DIR__, "../data/merge_20221012_local")

function local_opt(v, B, free_terms, free_term_limits)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    cs = VariableRef[]
    for (vf, (lb, ub)) in zip(free_terms, eachrow(free_term_limits))
        c = @variable(model)
        if isfinite(lb)
            @constraint(model, c >= lb)
        end
        if isfinite(ub)
            @constraint(model, c <= ub)
        end
        v = @expression(model, v .+ vf .* c)
        push!(cs, c)
    end
    nv, nt = size(B)
    @assert length(v) == nv
    t = @variable(model, [1:nt])
    v = @expression(model, B * t .+ v)
    maxv = @variable(model)
    @constraint(model, maxv .>= v)
    @constraint(model, maxv .>= .-v)
    @objective(model, Min, maxv)
    JuMP.optimize!(model)

    return value.(cs), value.(t), value.(v)
end

function solve_transfer(i)
    xpos_um = xpos_ums[i]
    @show xpos_um
    data = coeff_data[i]
    x0 = data["solution"]
    @assert data["xpos_um"] == xpos_um
    cs, ts, vs = local_opt(data["solution"],
                           data["free_solution"],
                           data["limited_solution"],
                           data["limits"])
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vs,
                "xpos_um"=>xpos_um,
                "limited_coeff"=>cs, "free_coeff"=>ts)
end
const transfer_solutions = [@time(solve_transfer(i)) for i in 1:length(xpos_ums)]

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
