#!/usr/bin/julia

include("gen_slice_data_utils.jl")

using NaCsPlot
using PyPlot

const loader = SliceLoader(ARGS[1])
slice_data = load(loader, ARGS[2], filter=get(ARGS, 3, nothing),
                  min_pos_um=-70, max_pos_um=70)

@show extrema(slice_data.voltages)

term_x1_raw = ideal_term(loader, slice_data.pos_um, (1, 0, 0),
                         min_pos_um=-70, max_pos_um=70)
term_x2_raw = ideal_term(loader, slice_data.pos_um, (2, 0, 0),
                         min_pos_um=-70, max_pos_um=70)
term_x3_raw = ideal_term(loader, slice_data.pos_um, (3, 0, 0),
                         min_pos_um=-70, max_pos_um=70)
term_x4_raw = ideal_term(loader, slice_data.pos_um, (4, 0, 0),
                         min_pos_um=-70, max_pos_um=70)

term_y1_raw = ideal_term(loader, slice_data.pos_um, (0, 1, 0),
                         min_pos_um=-70, max_pos_um=70)
term_y2_raw = ideal_term(loader, slice_data.pos_um, (0, 2, 0),
                         min_pos_um=-70, max_pos_um=70)

term_z1_raw = ideal_term(loader, slice_data.pos_um, (0, 0, 1),
                         min_pos_um=-70, max_pos_um=70)
term_z2_raw = ideal_term(loader, slice_data.pos_um, (0, 0, 2),
                         min_pos_um=-70, max_pos_um=70)

term_x1y1_raw = ideal_term(loader, slice_data.pos_um, (1, 1, 0),
                           min_pos_um=-70, max_pos_um=70)
term_y1z1_raw = ideal_term(loader, slice_data.pos_um, (0, 1, 1),
                           min_pos_um=-70, max_pos_um=70)
term_x1z1_raw = ideal_term(loader, slice_data.pos_um, (1, 0, 1),
                           min_pos_um=-70, max_pos_um=70)

term_x1y2_raw = ideal_term(loader, slice_data.pos_um, (1, 2, 0),
                           min_pos_um=-70, max_pos_um=70)
term_x2y2_raw = ideal_term(loader, slice_data.pos_um, (2, 2, 0),
                           min_pos_um=-70, max_pos_um=70)
term_x1z2_raw = ideal_term(loader, slice_data.pos_um, (1, 0, 2),
                           min_pos_um=-70, max_pos_um=70)
term_x2z2_raw = ideal_term(loader, slice_data.pos_um, (2, 0, 2),
                           min_pos_um=-70, max_pos_um=70)

term_dx = term_x1_raw * (Solutions.l_unit_um / Solutions.V_unit_uV)
term_dy = term_y1_raw * (Solutions.l_unit_um / Solutions.V_unit_uV)
term_dz = term_z1_raw * (Solutions.l_unit_um / Solutions.V_unit_uV)

term_xy = term_x1y1_raw
term_yz = term_y1z1_raw
term_zx = term_x1z1_raw
term_z2 = term_z2_raw - term_y2_raw
term_x2 = term_x2_raw - (term_y2_raw + term_z2_raw) / 2

term_x3 = term_x3_raw - (term_x1y2_raw + term_x1z2_raw) / 2
term_x4 = term_x4_raw - (term_x2y2_raw + term_x2z2_raw) / 2

target_term = nothing
if slice_data.name == "dx"
    target_term = term_dx
elseif slice_data.name == "dy"
    target_term = term_dy
elseif slice_data.name == "dz"
    target_term = term_dz
elseif slice_data.name == "xy"
    target_term = term_xy
elseif slice_data.name == "yz"
    target_term = term_yz
elseif slice_data.name == "zx"
    target_term = term_zx
elseif slice_data.name == "z2"
    target_term = term_z2
elseif slice_data.name == "x2"
    target_term = term_x2
elseif slice_data.name == "x3"
    target_term = term_x3
elseif slice_data.name == "x4"
    target_term = term_x4
end

figure(figsize=[6.4 * 3, 4.8 * 2])
subplot(2, 3, 1)
target_term !== nothing && plot(target_term.xpos, target_term.C, ls="--")
plot(slice_data.data.xpos, slice_data.data.C)
xlim([-70, 70])
grid()
title("\$1\$")

subplot(2, 3, 2)
target_term !== nothing && plot(target_term.xpos, target_term.Y, ls="--")
plot(slice_data.data.xpos, slice_data.data.Y)
xlim([-70, 70])
grid()
title("\$y\$")

subplot(2, 3, 3)
target_term !== nothing && plot(target_term.xpos, target_term.Z, ls="--")
plot(slice_data.data.xpos, slice_data.data.Z)
xlim([-70, 70])
grid()
title("\$z\$")

subplot(2, 3, 4)
target_term !== nothing && plot(target_term.xpos, target_term.Y2, ls="--")
plot(slice_data.data.xpos, slice_data.data.Y2)
xlim([-70, 70])
grid()
title("\$y^2\$")

subplot(2, 3, 5)
target_term !== nothing && plot(target_term.xpos, target_term.YZ, ls="--")
plot(slice_data.data.xpos, slice_data.data.YZ)
xlim([-70, 70])
grid()
title("\$yz\$")

subplot(2, 3, 6)
target_term !== nothing && plot(target_term.xpos, target_term.Z2, ls="--")
plot(slice_data.data.xpos, slice_data.data.Z2)
xlim([-70, 70])
grid()
title("\$z^2\$")

suptitle(slice_data.name)

NaCsPlot.maybe_show()
