#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using JuMP
using Ipopt
using MAT

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/merge_coeff_20221021.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const prefix = joinpath(@__DIR__, "../data/merge_20221021_local")

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

function solve_transfer(data)
    x0 = data["solution"]
    cs, ts, vs = local_opt(data["solution"],
                           data["free_solution"],
                           data["limited_solution"],
                           data["limits"])
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vs,
                "limited_coeff"=>cs, "free_coeff"=>ts)
end
const transfer_solutions = [[begin
                                 println("$i/$(length(section))")
                                 @time solve_transfer(data)
                             end for (i, data) in enumerate(section)]
                            for section in coeff_data]

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end
