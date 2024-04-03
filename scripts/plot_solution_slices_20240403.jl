#!/usr/bin/julia

include("gen_slice_data_utils.jl")

using NaCsPlot
using PyPlot

const loader = SliceDataLoader(ARGS[1])
slice_data = load(loader, ARGS[2], filter=get(ARGS, 3, nothing),
                  min_pos_um=-70, max_pos_um=70)

figure(figsize=[6.4 * 3, 4.8 * 2])
subplot(2, 3, 1)
plot(slice_data.data.xpos, slice_data.data.C)
xlim([-70, 70])
grid()
title("\$1\$")

subplot(2, 3, 2)
plot(slice_data.data.xpos, slice_data.data.Y)
xlim([-70, 70])
grid()
title("\$y\$")

subplot(2, 3, 3)
plot(slice_data.data.xpos, slice_data.data.Z)
xlim([-70, 70])
grid()
title("\$z\$")

subplot(2, 3, 4)
plot(slice_data.data.xpos, slice_data.data.Y2)
xlim([-70, 70])
grid()
title("\$y^2\$")

subplot(2, 3, 5)
plot(slice_data.data.xpos, slice_data.data.YZ)
xlim([-70, 70])
grid()
title("\$yz\$")

subplot(2, 3, 6)
plot(slice_data.data.xpos, slice_data.data.Z2)
xlim([-70, 70])
grid()
title("\$z^2\$")

suptitle(slice_data.name)

NaCsPlot.maybe_show()
