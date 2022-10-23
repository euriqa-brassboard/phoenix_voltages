#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

import PhoenixVoltages.Optimizers
using MAT
import Ipopt
using JuMP

const coeff_data, electrode_names = matopen(joinpath(@__DIR__, "../data/merge_coeff_20221021.mat")) do mat
    return read(mat, "data"), read(mat, "electrode_names")
end

const prev_solution = matopen(ARGS[1]) do mat
    @assert read(mat, "electrode_names") == electrode_names
    return read(mat, "transfer_solutions")
end

const nsection = length(coeff_data)
const n_sol_points = length.(coeff_data)

@assert all(length.(prev_solution) == n_sol_points)

struct SolutionInfo
    maxv::Float64
    maxvdiff::Float64
end
const prev_solution_cache = Vector{SolutionInfo}.(undef, n_sol_points)
function fill_solutions_cache(si)
    local prev_vmap
    for i in 1:n_sol_points[si]
        data = coeff_data[si][i]
        sol = prev_solution[si][i]
        electrodes = data["electrodes"]
        @assert sol["electrodes"] == electrodes
        v = data["solution"] .+ data["free_solution"] * sol["free_coeff"]
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
        prev_solution_cache[si][i] = SolutionInfo(maxv, maxvdiff)
    end
end
for si in 1:nsection
    fill_solutions_cache(si)
end

mutable struct ModelSection
    cs::Vector{Vector{VariableRef}}
    ts::Vector{Vector{VariableRef}}

    # Electrode voltages
    vs::Vector{Vector{AffExpr}}

    # Maximum voltage
    maxvs::Vector{VariableRef}
    # Maximum voltage difference with last frame
    maxvdiffs::Vector{VariableRef}

    # Block cost functions
    maxmaxvs::Vector{VariableRef}
    maxmaxvdiffs::Vector{VariableRef}

    first_vmap::Dict{Int,AffExpr}
    last_vmap::Dict{Int,AffExpr}

    function ModelSection()
        return new(Vector{VariableRef}[], Vector{VariableRef}[],
                   Vector{AffExpr}[], VariableRef[], VariableRef[],
                   VariableRef[], VariableRef[])
    end
end

struct ModelBuilder
    model::Model
    sections::Vector{ModelSection}

    function ModelBuilder()
        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "max_iter",
                                parse(Int, get(ENV, "OPT_MAX_ITER", "30000")))
        set_optimizer_attribute(model, "print_level", 5)
        return new(model, [ModelSection() for i in 1:nsection])
    end
end

struct TrapWeights
    maxvs::Float64
    maxvdiffs::Float64
    maxmaxv::Float64
    maxmaxvdiff::Float64
    hard_maxv::Float64
end

function maxdiff_var(model, vmap1, vmap2)
    maxvdiff = @variable(model)
    for k in union(keys(vmap1), keys(vmap2))
        v1 = get(vmap1, k, 0)
        v2 = get(vmap2, k, 0)
        @constraint(model, maxvdiff >= v1 - v2)
        @constraint(model, maxvdiff >= v2 - v1)
    end
    return maxvdiff
end

function add_block!(builder::ModelBuilder, si, from, to)
    model = builder.model
    section = builder.sections[si]
    maxmaxv = @variable(model)
    maxmaxv_init = 0.0
    maxmaxvdiff = @variable(model)
    maxmaxvdiff_init = 0.0
    for i in from:to
        println("$i/$(length(coeff_data[si]))")
        data = coeff_data[si][i]
        sol = prev_solution[si][i]
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
        push!(section.cs, cs)
        B = data["free_solution"]
        t_init = sol["free_coeff"]
        nv, nt = size(B)
        @assert length(v) == nv
        @assert length(t_init) == nt

        t = @variable(model, [j=1:nt], start=t_init[j])
        maxv_init = prev_solution_cache[si][i].maxv
        maxv = @variable(model, start=maxv_init)
        push!(section.ts, t)

        v = @expression(model, B * t .+ v)
        push!(section.vs, v)
        @constraint(model, maxv .>= v)
        @constraint(model, maxv .>= .-v)
        push!(section.maxvs, maxv)

        vmap = Dict(zip(electrodes, v))
        if !isdefined(section, :first_vmap)
            section.first_vmap = vmap
        end
        if isdefined(section, :last_vmap)
            maxvdiff = maxdiff_var(model, vmap, section.last_vmap)
            maxvdiff_init = prev_solution_cache[si][i].maxvdiff
            set_start_value(maxvdiff, maxvdiff_init)
            push!(section.maxvdiffs, maxvdiff)
            @constraint(model, maxmaxvdiff >= maxvdiff)
            maxmaxvdiff_init = max(maxmaxvdiff_init, maxvdiff_init)
        end
        section.last_vmap = vmap

        @constraint(model, maxmaxv >= maxv)
        maxmaxv_init = max(maxmaxv_init, maxv_init)
    end
    set_start_value(maxmaxv, maxmaxv_init)
    set_start_value(maxmaxvdiff, maxmaxvdiff_init)

    push!(section.maxmaxvs, maxmaxv)
    push!(section.maxmaxvdiffs, maxmaxvdiff)

    return builder
end

function add_section!(builder::ModelBuilder, si, block_sz)
    for i in 1:block_sz:n_sol_points[si]
        add_block!(builder, si, i, min(i + block_sz - 1, n_sol_points[si]))
    end
end

function add_all_sections!(builder::ModelBuilder, block_sz)
    for si in 1:nsection
        add_section!(builder, si, block_sz)
    end
    ### Hard code inter-section diff relation
    model = builder.model
    sections = builder.sections

    merge_enter_diff = maxdiff_var(model, sections[1].last_vmap, sections[2].first_vmap)
    move_enter_diff = maxdiff_var(model, sections[1].last_vmap, sections[3].first_vmap)
    merge_end_diff = maxdiff_var(model, sections[2].last_vmap, sections[3].last_vmap)

    push!(sections[1].maxvdiffs, merge_enter_diff)
    push!(sections[1].maxvdiffs, move_enter_diff)
    push!(sections[2].maxvdiffs, merge_end_diff)

    @constraint(model, sections[1].maxmaxvdiffs[end] >= merge_enter_diff)
    @constraint(model, sections[1].maxmaxvdiffs[end] >= move_enter_diff)
    @constraint(model, sections[2].maxmaxvdiffs[1] >= merge_enter_diff)
    @constraint(model, sections[3].maxmaxvdiffs[1] >= move_enter_diff)
    @constraint(model, sections[2].maxmaxvdiffs[end] >= merge_end_diff)
    @constraint(model, sections[3].maxmaxvdiffs[end] >= merge_end_diff)
    return
end

function section_obj(builder::ModelBuilder, section::ModelSection,
                     block_sz, weights::TrapWeights)
    model = builder.model
    obj = sum(section.maxvs) * weights.maxvs
    if weights.maxvdiffs > 0
        obj += sum(section.maxvdiffs) * weights.maxvdiffs
    end
    if weights.maxmaxv > 0
        obj += sum(section.maxmaxvs) * (weights.maxmaxv * block_sz)
    end
    if weights.maxmaxvdiff > 0
        obj += sum(section.maxmaxvdiffs) * (weights.maxmaxvdiff * block_sz)
    end
    if isfinite(weights.hard_maxv)
        for v in section.maxmaxvs
            @constraint(model, v <= weights.hard_maxv)
        end
    end
    return obj
end

function finalize_model!(builder::ModelBuilder, block_sz, weights::TrapWeights)
    model = builder.model
    obj = sum(section_obj(builder, section, block_sz, weights)
              for section in builder.sections)
    @objective(model, Min, obj)
    return builder
end

function optimize_model!(builder::ModelBuilder)
    @time JuMP.optimize!(builder.model)
    return [([value.(v) for v in section.cs], [value.(v) for v in section.ts],
             [value.(v) for v in section.vs]) for section in builder.sections]
end

const prefix = ARGS[2]
const block_size = parse(Int, ARGS[3])

function pack_data(data, cs, ts, vs)
    return Dict("electrodes"=>data["electrodes"],
                "voltages"=>vs, "limited_coeff"=>cs, "free_coeff"=>ts)
end

function run()
    builder = ModelBuilder()
    add_all_sections!(builder, block_size)
    finalize_model!(builder, block_size, TrapWeights(0.2, 0.1, 0.1, 1.5, 29))
    res = optimize_model!(builder)
    transfer_solutions = [[pack_data(data, c, t, v) for (data, c, t, v)
                               in zip(cd, cs, ts, vs)]
                          for (cd, (cs, ts, vs)) in zip(coeff_data, res)]
    matopen("$(prefix).mat", "w") do mat
        write(mat, "electrode_names", electrode_names)
        write(mat, "transfer_solutions", transfer_solutions)
    end
end

run()
