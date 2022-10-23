#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using Printf
using PhoenixVisual
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions

function load_solution(name)
    return matopen(name) do mat
        return (transfer_solutions=read(mat, "transfer_solutions"),
                electrode_names=read(mat, "electrode_names"))
    end
end

const transfer = load_solution(ARGS[2])
const electrode_names = Vector{String}.(transfer.electrode_names)

const coeff_data = matopen(ARGS[3]) do mat
    return read(mat, "data")
end

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file,
                                              electrode_names=electrode_names)
const fits_cache = Potentials.FitCache(Fitting.PolyFitter(2, 2, 2), solution)
const centers = Solutions.CenterTracker()

function get_xpos_um(data)
    xpos_um_min = [minimum(d["xpos_um"]) for d in data]
    xpos_um_max = [maximum(d["xpos_um"]) for d in data]
    xpos_min_diff = abs(xpos_um_min[1] - xpos_um_min[end])
    xpos_max_diff = abs(xpos_um_max[1] - xpos_um_max[end])
    return xpos_max_diff > xpos_min_diff ? xpos_um_max : xpos_um_min
end
const xpos_ums = get_xpos_um.(coeff_data)

const prefix = ARGS[4]

const xrange_ums = -3220:1330

function create_frame(title, sol, data)
    xpos_ums = data["xpos_um"]
    frame = PhoenixVisual.render_frame(fits_cache, centers,
                                       sol["electrodes"], sol["voltages"] ./ 1.5,
                                       title=title,
                                       xpos_um=xpos_ums[1], plotx_ums=xrange_ums,
                                       finalize=false)
    for i in 2:length(xpos_ums)
        xpos_um = xpos_ums[i]
        xidx = Potentials.x_axis_to_index(solution, xpos_um)
        yidx = get(centers, xidx)[1]
        ypos_um = Potentials.y_index_to_axis(solution, yidx) * 1000
        c = PhoenixVisual.add_circle!(frame, xpos_um, ypos_um, 20)
        c["fill"] = "blue"
    end
    PhoenixVisual.finalize_svg!(frame)
    return frame
end

const xpos_prefixes = ["Single ion", "Chain", "Chain"]

for si in 1:length(coeff_data)
    solutions = transfer.transfer_solutions[si]
    cd = coeff_data[si]
    xpos_um = xpos_ums[si]
    xpos_prefix = xpos_prefixes[si]
    mkpath("$(prefix)/$si")
    for (i, (solution, data)) in enumerate(zip(solutions, cd))
        println("$i/$(length(cd))")
        template = create_frame(@sprintf("%s: % 5d μm", xpos_prefix, xpos_um[i]),
                                solution, data)
        open("$(prefix)/$si/$i.svg", "w") do io
            write(io, template)
        end
    end
end
