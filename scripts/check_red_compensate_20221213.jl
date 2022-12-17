#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Reconstruct
using PhoenixVoltages.Outputs
using MAT
using PhoenixVisual

const potential_file = ARGS[1]
const potential = Potentials.import_pillbox_64(potential_file)
const fits_cache = Solutions.compensate_fitter1(potential, sizes=(5, 5, 15))
const centers = Solutions.CenterTracker()

function get_all_fits(xpos_um, voltages, electrode_names)
    electrodes_voltages = Dict{String,Float64}(zip(electrode_names, voltages))
    xidx = Potentials.x_axis_to_index(potential, xpos_um ./ 1000)
    cf = Reconstruct.get_center_fit(fits_cache, electrodes_voltages,
                                    xidx, centers=centers)
    return Solutions.get_compensate_terms1(cf, potential.stride .* 1000)
end

function fix_electrode_name(orig_name)
    m = match(r"^Q0([0-9])$", orig_name)
    if m !== nothing
        return "Q$(m[1])"
    end
    m = match(r"^S0([0-9])$", orig_name)
    if m !== nothing
        return "S$(m[1])"
    end
    return orig_name
end

function show_file(file::CompensationFile)
    file.map.names .= fix_electrode_name.(file.map.names)
    for (name, values) in zip(file.term_names, file.term_values)
        if startswith(name, "Q")
            continue
        end
        println("$(strip(name)):")
        println("    ", get_all_fits(-175.0, values, file.map.names))
    end
end

function show_file(file::TransferFile)
    file.map.names .= fix_electrode_name.(file.map.names)
    println("Trap solution:")
    println("    ", get_all_fits(-175.0, file.line_values[1], file.map.names))
    return
end

try
    show_file(load_file(ARGS[2], CompensationFile))
catch
    show_file(load_file(ARGS[2], TransferFile))
end
