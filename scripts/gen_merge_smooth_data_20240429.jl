#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT
import Ipopt
using JuMP

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
    load_coeff_multi(joinpath(@__DIR__, "../data/merge_coeff_20240429"))

const prev_solution = matopen(ARGS[1]) do mat
    @assert read(mat, "electrode_names") == electrode_names
    return read(mat, "transfer_solutions")
end

const n_sol_points = length(coeff_data)

const xpos_ums = [data["xpos_um"] for data in coeff_data]

@assert length(prev_solution) == n_sol_points
@assert [data["xpos_um"] for data in prev_solution] == xpos_ums

struct SolutionInfo
    maxv::Float64
    maxvdiff::Float64
end
const prev_solution_cache = Vector{SolutionInfo}(undef, n_sol_points)
function fill_solutions_cache()
    local prev_vmap
    for i in 1:n_sol_points
        data = coeff_data[i]
        sol = prev_solution[i]
        electrodes = data["electrodes"]
        @assert sol["electrodes"] == electrodes
        v = copy(data["solution"])
        for (vf, c) in zip(data["limited_solution"], sol["limited_coeff"])
            v .= v .+ vf .* c
        end
        vmap = Dict(zip(electrodes, v))
        maxvdiff = 0.0
        if i > 1
            for k in union(keys(vmap), keys(prev_vmap))
                v1 = get(vmap, k, 0)
                v2 = get(prev_vmap, k, 0)
                maxvdiff = max(maxvdiff, abs(v1 - v2))
            end
        end
        prev_vmap = vmap
        maxv = maximum(abs.(v))
        prev_solution_cache[i] = SolutionInfo(maxv, maxvdiff)
    end
end
fill_solutions_cache()

mutable struct ModelBuilder
    model::Model
    cs::Vector{Vector{VariableRef}}

    # Electrode voltages
    vs::Vector{Vector{AffExpr}}

    # Maximum voltage
    maxvs::Vector{VariableRef}
    # Maximum voltage difference with last frame
    maxvdiffs::Vector{VariableRef}

    # Block cost functions
    maxmaxvs::Vector{VariableRef}
    maxmaxvdiffs::Vector{VariableRef}

    prev_vmap::Dict{Int,AffExpr}
    function ModelBuilder()
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "max_iter",
                                parse(Int, get(ENV, "OPT_MAX_ITER", "30000")))
        set_optimizer_attribute(model, "print_level", 5)
        return new(model, Vector{VariableRef}[],
                   Vector{AffExpr}[], VariableRef[], VariableRef[],
                   VariableRef[], VariableRef[])
    end
end

struct TrapWeights
    maxvs::Float64
    maxvdiffs::Float64
    maxmaxv::Float64
    maxmaxvdiff::Float64
    hard_maxv::Float64
end

function add_block!(builder::ModelBuilder, from, to)
    model = builder.model
    maxmaxv = @variable(model)
    maxmaxv_init = 0.0
    maxmaxvdiff = @variable(model)
    maxmaxvdiff_init = 0.0
    for i in from:to
        data = coeff_data[i]
        # @show data["xpos_um"]
        sol = prev_solution[i]
        electrodes = data["electrodes"]
        v = data["solution"]
        cs = VariableRef[]
        for (vf, (lb, ub), ci) in zip(data["limited_solution"], eachrow(data["limits"]),
                                      sol["limited_coeff"])
            c = @variable(model, start=ci)
            push!(cs, c)
            if isfinite(lb)
                @constraint(model, c >= lb)
            end
            if isfinite(ub)
                @constraint(model, c <= ub)
            end
            v = @expression(model, v .+ vf .* c)
        end
        push!(builder.cs, cs)

        maxv_init = prev_solution_cache[i].maxv
        maxv = @variable(model, start=maxv_init)
        push!(builder.vs, v)


        for (cterm, (lb, ub)) in zip(data["constraints_terms"],
                                     eachrow(data["constraints_limits"]))
            cterm_val = @expression(model, cterm[1] + sum(cs .* cterm[2:end]))
            if isfinite(lb)
                @constraint(model, cterm_val >= lb)
            end
            if isfinite(ub)
                @constraint(model, cterm_val <= ub)
            end
        end

        @constraint(model, maxv .>= v)
        @constraint(model, maxv .>= .-v)
        push!(builder.maxvs, maxv)

        vmap = Dict(zip(electrodes, v))
        if isdefined(builder, :prev_vmap)
            maxvdiff_init = prev_solution_cache[i].maxvdiff
            maxvdiff = @variable(model, start=maxvdiff_init)
            prev_vmap = builder.prev_vmap
            for k in union(keys(vmap), keys(prev_vmap))
                v1 = get(vmap, k, 0)
                v2 = get(prev_vmap, k, 0)
                @constraint(model, maxvdiff >= v1 - v2)
                @constraint(model, maxvdiff >= v2 - v1)
            end
            push!(builder.maxvdiffs, maxvdiff)
            @constraint(model, maxmaxvdiff >= maxvdiff)
            maxmaxvdiff_init = max(maxmaxvdiff_init, maxvdiff_init)
        end
        builder.prev_vmap = vmap

        @constraint(model, maxmaxv >= maxv)
        maxmaxv_init = max(maxmaxv_init, maxv_init)
    end
    set_start_value(maxmaxv, maxmaxv_init)
    set_start_value(maxmaxvdiff, maxmaxvdiff_init)

    push!(builder.maxmaxvs, maxmaxv)
    push!(builder.maxmaxvdiffs, maxmaxvdiff)

    return builder
end

function add_all_blocks!(builder::ModelBuilder, block_sz)
    for i in 1:block_sz:n_sol_points
        add_block!(builder, i, min(i + block_sz - 1, n_sol_points))
    end
end

function finalize_model!(builder::ModelBuilder, block_sz, weights::TrapWeights)
    model = builder.model

    obj = sum(builder.maxvs) * weights.maxvs
    if weights.maxvdiffs > 0
        obj += sum(builder.maxvdiffs) * weights.maxvdiffs
    end
    if weights.maxmaxv > 0
        obj += sum(builder.maxmaxvs) * (weights.maxmaxv * block_sz)
    end
    if weights.maxmaxvdiff > 0
        obj += sum(builder.maxmaxvdiffs) * (weights.maxmaxvdiff * block_sz)
    end
    if isfinite(weights.hard_maxv)
        for v in builder.maxmaxvs
            @constraint(model, v <= weights.hard_maxv)
        end
    end
    @objective(model, Min, obj)
    return builder
end

function optimize_model!(builder::ModelBuilder)
    @time JuMP.optimize!(builder.model)
    return ([value.(v) for v in builder.cs], [value.(v) for v in builder.vs])
end

const prefix = ARGS[2]
const block_size = parse(Int, ARGS[3])

function pack_data(data, cs, vs)
    return Dict("electrodes"=>data["electrodes"],
                "voltages"=>vs, "limited_coeff"=>cs, "xpos_um"=>data["xpos_um"])
end

function run()
    builder = ModelBuilder()
    add_all_blocks!(builder, block_size)
    finalize_model!(builder, block_size, TrapWeights(0.2, 0.1, 0.1, 1.5, 199))
    cs, vs = optimize_model!(builder)
    transfer_solutions = [pack_data(data, c, v) for (data, c, v)
                              in zip(coeff_data, cs, vs)]
    matopen("$(prefix).mat", "w") do mat
        write(mat, "electrode_names", electrode_names)
        write(mat, "transfer_solutions", transfer_solutions)
    end
end

run()
