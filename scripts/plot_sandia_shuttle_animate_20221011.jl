#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using MAT
using Printf
using PhoenixVisual
using PhoenixVoltages.Outputs
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Solutions

const solution_file = ARGS[1]
const solution = Potentials.import_pillbox_64(solution_file)
const fits_cache = Potentials.FitCache(Fitting.PolyFitter(2, 2, 2), solution)

const transfer_file = load_file(ARGS[2], TransferFile)
const electrodes_idx = [get(solution.electrode_index, name, 1)
                        for name in transfer_file.map.names]

const prefix = ARGS[3]

const xrange_ums = -3220:1330

function create_frame(i, data)
    xpos_um = i - 176 - 3045
    @show i, xpos_um
    return PhoenixVisual.render_frame(fits_cache, Solutions.CenterTracker(),
                                      electrodes_idx, data .* 100,
                                      title=@sprintf("% 5d μm", xpos_um),
                                      xpos_um=xpos_um, plotx_ums=xrange_ums)
end

mkpath("$(prefix)")
for (i, data) in enumerate(transfer_file.line_values[1:end - 1])
    template = create_frame(i, data)
    open("$(prefix)/$i.svg", "w") do io
        write(io, template)
    end
end
