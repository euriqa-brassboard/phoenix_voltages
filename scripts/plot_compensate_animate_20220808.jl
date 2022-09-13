#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using Printf
using Statistics
using PhoenixVisual
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions

function load_solution(name)
    return matopen(name) do mat
        return (transfer_solutions=read(mat, "transfer_solutions"),
                electrode_names=read(mat, "electrode_names"),
                termname=read(mat, "termname"))
    end
end

const centers = matopen(joinpath(@__DIR__, "../data/rf_center.mat")) do mat
    return Solutions.CenterTracker(read(mat, "zy_index"))
end
const transfer = load_solution(ARGS[2])
const electrode_names = Vector{String}.(transfer.electrode_names)

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file,
                                              electrode_names=electrode_names)
const fits_cache = Potentials.FitCache(Fitting.PolyFitter(2, 2, 2), solution)

const prefix = ARGS[3]

const xrange_ums = -3220:1330

function get_scaling()
    m = median(maximum(abs.(data["voltages"]))
               for data in transfer.transfer_solutions if -500 < data["xpos_um"] < 500)
    return 10 / m
end

function get_offset()
    if transfer.termname in ("dx", "x3", "xy", "zx")
        return 0.5
    else
        return 0.2
    end
end

const yscale = get_scaling()
const yoffset = get_offset()

function create_frame(data)
    xpos_um = data["xpos_um"]
    title = transfer.termname * @sprintf(" @ % 5d μm", xpos_um)
    return PhoenixVisual.render_frame(fits_cache, centers,
                                      data["electrodes"], data["voltages"],
                                      title=title, xpos_um=xpos_um,
                                      plotx_ums=xrange_ums, plot_yoffset=yoffset)
end

mkpath("$(prefix)_$(transfer.termname)")
for (i, data) in enumerate(transfer.transfer_solutions)
    @show i data["xpos_um"]
    template = create_frame(data)
    open("$(prefix)_$(transfer.termname)/$i.svg", "w") do io
        write(io, template)
    end
end
