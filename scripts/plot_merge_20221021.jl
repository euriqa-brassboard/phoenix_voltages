#!/usr/bin/julia

using MAT
using PyPlot
using NaCsPlot

term_to_map(t) = Dict(name=>value for (name, value) in zip(t["electrodes"], t["voltages"]))
function vdiff(t1, t2)
    map1 = term_to_map(t1)
    map2 = term_to_map(t2)
    maxdiff = 0.0
    for k in union(keys(map1), keys(map2))
        v1 = get(map1, k, 0.0)
        v2 = get(map2, k, 0.0)
        d = abs(v1 - v2)
        maxdiff = max(d, maxdiff)
    end
    return maxdiff
end

const coeff_data = matopen(joinpath(@__DIR__, "../data/merge_coeff_20221021.mat")) do mat
    return read(mat, "data")
end
const nsection = length(coeff_data)
const plot_ncol = ceil(Int, sqrt(nsection))
const plot_nrow = ceil(Int, nsection / plot_ncol)

function get_xpos_um(data)
    xpos_um_min = [minimum(d["xpos_um"]) for d in data]
    xpos_um_max = [maximum(d["xpos_um"]) for d in data]
    xpos_min_diff = abs(xpos_um_min[1] - xpos_um_min[end])
    xpos_max_diff = abs(xpos_um_max[1] - xpos_um_max[end])
    return xpos_max_diff > xpos_min_diff ? xpos_um_max : xpos_um_min
end
const xpos_ums = get_xpos_um.(coeff_data)

function load_section(si, section)
    maxv = [maximum(abs.(sol["voltages"])) for sol in section]
    maxdiff = [vdiff(section[i + 1], section[i]) for i in 1:length(section) - 1]
    return (xpos_um=xpos_ums[si], maxv=maxv, maxdiff=maxdiff)
end

function load_solution(name)
    solution = matopen(joinpath(@__DIR__, "../data", name)) do mat
        return read(mat, "transfer_solutions")
    end
    @assert length(solution) == nsection
    return [load_section(i, section) for (i, section) in enumerate(solution)]
end

const solutions = [("local", load_solution("merge_20221021_local.mat")),
                   # ("blk 70", load_solution("merge_20221021_blk70.mat"))
                   ]

const prefix = joinpath(@__DIR__, "../imgs/merge_20221021")

figure(figsize=[6.4 * plot_ncol, 4.8 * plot_nrow])
for si in 1:nsection
    subplot(plot_nrow, plot_ncol, si)
    for (label, solution) in solutions
        plot(solution[si].xpos_um, solution[si].maxv, label=label)
    end
    # ylim([0.5, 30])
    legend(fontsize=10, ncol=3)
    grid()
    xlabel("Position (μm)")
    ylabel("Max (V)")
end
NaCsPlot.maybe_save("$(prefix)_max")

figure(figsize=[6.4 * plot_ncol, 4.8 * plot_nrow])
for si in 1:nsection
    subplot(plot_nrow, plot_ncol, si)
    for (label, solution) in solutions
        plot(solution[si].xpos_um[1:end - 1], solution[si].maxdiff, label=label)
    end
    # ylim([0, 2.5])
    legend(fontsize=10, ncol=3)
    grid()
    xlabel("Position (μm)")
    ylabel("Difference (V)")
end
NaCsPlot.maybe_save("$(prefix)_diff")

NaCsPlot.maybe_show()
