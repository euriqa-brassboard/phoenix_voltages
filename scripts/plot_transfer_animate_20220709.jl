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

function get_rf_center(xpos_um)
    xidx = Potentials.x_axis_to_index(solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

const xrange_ums = -3220:1330

function create_frame(data)
    template = PhoenixVisual.get_template(plot=true)

    xpos_um = data["xpos_um"]
    rf_center_idx = get_rf_center(xpos_um)
    ypos_um = Potentials.y_index_to_axis(solution, rf_center_idx[2]) * 1000

    PhoenixVisual.set_title!(template, @sprintf("% 5d μm", xpos_um))

    c = PhoenixVisual.add_circle!(template, xpos_um, ypos_um, 20)
    c["fill"] = "blueviolet"

    electrode_voltages = zip(data["electrodes"], data["voltages"])
    function get_voltage(x_um)
        local pos = get_rf_center(x_um)
        fit = Potentials.get_multi_electrodes(fits_cache, electrode_voltages,
                                              (pos[3], pos[2], pos[1]))
        return fit[0, 0, 0]
    end
    ax_potential = get_voltage.(xrange_ums)

    plot_yoffset = 0.2
    line2 = PhoenixVisual.add_plotline!(template,
                                        [xrange_ums[1] - 200, xrange_ums[end] + 200],
                                        [plot_yoffset, plot_yoffset])
    line2["fill"] = "none"
    line2["stroke"] = "gainsboro"
    line = PhoenixVisual.add_plotline!(template, xrange_ums,
                                       ax_potential .* 0.2 .+ plot_yoffset)
    line["fill"] = "none"
    line["stroke"] = "coral"

    voltage_map = Dict{String,Float64}()
    for (idx, v) in electrode_voltages
        v = v / 20
        for name in electrode_names[idx]
            voltage_map[name] = v
        end
    end

    PhoenixVisual.fill_electrodes!(template, voltage_map)

    PhoenixVisual.finalize_svg!(template)
    return template
end

mkpath("$(prefix)")
for (i, data) in enumerate(transfer.transfer_solutions)
    @show i data["xpos_um"]
    template = create_frame(data)
    open("$(prefix)/$i.svg", "w") do io
        write(io, template)
    end
end
