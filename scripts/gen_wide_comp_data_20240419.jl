#!/usr/bin/julia

using MAT

using JuMP
using Ipopt

const electrode_names = matopen(joinpath(@__DIR__, "../data/wide_comp_coeff_20240419.mat")) do mat
    return read(mat, "electrode_names")
end

const xpos_ums = -350:350
const coeff_data = [matopen(joinpath(@__DIR__, "../data/wide_comp_coeff_20240419/$(xpos_um).mat")) do mat
                        return read(mat, "data")
                    end for xpos_um in xpos_ums]

const outputdir = joinpath(@__DIR__, "../data/wide_comp_20240419")

function fit_term(model::Model, coeff_data, term_idx, maxv, yz_weight;
                  relax_high_order=false)
    center_x0 = coeff_data["solution"]

    x0 = @view(center_x0[:, term_idx])
    B = coeff_data["free_solution"]
    nx, nt = size(B)
    if true
        @variable(model, t[1:nt])
        x = @expression(model, B * t .+ x0)
        if relax_high_order
            nhigh_order = size(center_x0, 2) - 10
            @variable(model, ht[1:nhigh_order])
            x = @expression(model, x .+ center_x0[:, 11:end] * ht)
        end
        @constraint(model, -maxv .<= x .<= maxv)
    else
        x = x0
    end

    nvars = nx + 1
    @variable(model, dc_offset)
    vars = [dc_offset; x]
    target = @view(coeff_data["slice_terms"][:, term_idx])
    err = @expression(model, coeff_data["slice_coeff"] * vars .- target)
    nvals = length(err)
    @variable(model, abserr[1:nvals])
    @constraint(model, .-abserr .<= err)
    @constraint(model, abserr .>= err)
    npoints = length(target) ÷ 6
    abserr = reshape(abserr, (6, npoints))

    @variable(model, maxerr)

    point_center = (1 + npoints) / 2
    max_offset = point_center - 1
    C = log(2) / max_offset^2
    for i in 1:npoints
        x_weight = exp((i - point_center)^2 * C)
        for j in 1:6
            e = abserr[j, i] * x_weight * yz_weight[j]
            @constraint(model, -maxerr <= e)
            @constraint(model, maxerr >= e)
        end
    end
    @objective(model, Min, maxerr)
    optimize!(model)
    return value.(vars)
end

function solve_term(coeff_data, term_idx)
    @show coeff_data["xpos_um"]
    model = Model(Ipopt.Optimizer)
    set_attribute(model, "print_level", 0)
    relax_high_order=false
    if term_idx == 1
        # dx
        maxv = 0.003
        yz_weight = (10, 3, 3, 1, 3, 1)
    elseif term_idx == 2
        # dy
        maxv = 0.003
        yz_weight = (5, 6, 3, 1, 3, 1)
    elseif term_idx == 3
        # dz
        maxv = 0.003
        yz_weight = (5, 5, 10, 1, 3, 1)
    elseif term_idx == 4
        # xy
        maxv = 6.0
        yz_weight = (400, 300, 6, 20, 200, 20)
    elseif term_idx == 5
        maxv = 3.2
        yz_weight = (40, 30, 60, 20, 300, 20)
        relax_high_order = true
    elseif term_idx == 6
        # zx
        maxv = 9.0
        yz_weight = (4, 30, 60, 40, 60, 40)
        relax_high_order = true
    elseif term_idx == 7
        # z2
        maxv = 10.0
        yz_weight = (40, 6, 6, 200, 300, 2000)
    elseif term_idx == 8
        # x2
        maxv = 19.0
        yz_weight = (40, 30, 3, 300, 300, 2000)
    elseif term_idx == 9
        # x3
        maxv = 200.0
        yz_weight = (6, 150, 20, 300, 700, 200)
        relax_high_order = true
    elseif term_idx == 10
        # x4
        maxv = 6000.0
        yz_weight = (4, 10, 10, 300, 300, 200)
        relax_high_order = true
    else
        error("Unknown term index: $term_idx")
    end
    return fit_term(model, coeff_data, term_idx, maxv, yz_weight;
                    relax_high_order=relax_high_order)
end

const term_names = ("dx", "dy", "dz", "xy", "yz", "zx",
                    "z2", "x2", "x3", "x4")

function pack_data(data, vals)
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end

function solve_all(term_idx)
    return solve_term.(coeff_data, term_idx)
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
