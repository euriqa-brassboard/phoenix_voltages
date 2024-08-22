#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using JuMP
using Ipopt
using MAT

function load_coeff_multi(dir)
    datas = Tuple{Int,Vector{Any}}[]
    for name in readdir(dir)
        m = match(r"data(\d+)\.mat", name)
        m === nothing && continue
        push!(datas, (parse(Int, m[1]), matread(joinpath(dir, name))["data"]))
    end
    sort!(datas, by=x->x[1])
    return (reduce(vcat, (data[2] for data in datas)),
            matread(joinpath(dir, "electrode_names.mat"))["electrode_names"])
end

const coeff_data, electrode_names =
    load_coeff_multi(joinpath(@__DIR__, "../data/merge_coeff_20240815"))

const xpos_ums = [data["xpos_um"] for data in coeff_data]

const prefix = joinpath(@__DIR__, "../data/merge_20240815_local")

function local_opt(v, free_terms, free_term_limits,
                   constraints_terms, constraints_limits)
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
    maxv = @variable(model)
    @constraint(model, maxv .>= v)
    @constraint(model, maxv .>= .-v)

    constraints = Tuple{AffExpr,Float64,Float64}[]
    for (cterm, (lb, ub)) in zip(constraints_terms, eachrow(constraints_limits))
        cterm_val = @expression(model, cterm[1] + sum(cs .* cterm[2:end]))
        push!(constraints, (cterm_val, lb, ub))
        if isfinite(lb)
            @constraint(model, cterm_val >= lb)
        end
        if isfinite(ub)
            @constraint(model, cterm_val <= ub)
        end
    end

    @objective(model, Min, maxv)
    JuMP.optimize!(model)

    for (cterm_val, lb, ub) in constraints
        cterm_val = value(cterm_val)
        if cterm_val > ub + 2e-5 || cterm_val < lb - 2e-5
            @warn "Constraint value $cterm_val outside of range [$lb, $ub]"
        end
    end

    @show value(maxv)
    if value(maxv) > 199
        @warn "Maximum voltage $(value(maxv)) too high."
    end

    return value.(cs), value.(v)
end

function solve_transfer(i)
    xpos_um = xpos_ums[i]
    @show xpos_um
    data = coeff_data[i]
    x0 = data["solution"]
    @assert data["xpos_um"] == xpos_um
    cs, vs = local_opt(data["solution"], data["limited_solution"], data["limits"],
                       data["constraints_terms"], data["constraints_limits"])
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vs,
                "xpos_um"=>xpos_um, "limited_coeff"=>cs)
end
const transfer_solutions = [@time(solve_transfer(i)) for i in 1:length(xpos_ums)]

matopen("$(prefix).mat", "w") do mat
    write(mat, "electrode_names", electrode_names)
    write(mat, "transfer_solutions", transfer_solutions)
end