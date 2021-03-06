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

function load_solution(name)
    solution = matopen(joinpath(@__DIR__, "../data", name)) do mat
        return read(mat, "transfer_solutions")
    end
    xpos_um = [sol["xpos_um"] for sol in solution]
    maxv = [maximum(abs.(sol["voltages"])) for sol in solution]
    maxdiff = [vdiff(solution[i + 1], solution[i]) for i in 1:length(solution) - 1]
    return (xpos_um=xpos_um, maxv=maxv, maxdiff=maxdiff)
end

const solution_ind = load_solution("transfer_20220705.mat")
const solution_smooth = load_solution("transfer_smooth_20220705.mat")
const solution_smooth_noglobal = load_solution("transfer_smooth_20220705_noglobal.mat")

const prefix = joinpath(@__DIR__, "../imgs/transfer_20220705")

figure()
plot(solution_ind.xpos_um, solution_ind.maxv, label="Individual")
plot(solution_smooth_noglobal.xpos_um, solution_smooth_noglobal.maxv, label="Neighbor")
plot(solution_smooth.xpos_um, solution_smooth.maxv, label="Global")
ylim([8, 17])
legend(fontsize=10, ncol=3)
grid()
xlabel("Position (μm)")
xlabel("Max (V)")
NaCsPlot.maybe_save("$(prefix)_max")

figure()
plot(solution_ind.xpos_um[1:end - 1], solution_ind.maxdiff, label="Individual")
plot(solution_smooth_noglobal.xpos_um[1:end - 1], solution_smooth_noglobal.maxdiff,
     label="Neighbor")
plot(solution_smooth.xpos_um[1:end - 1], solution_smooth.maxdiff, label="Global")
ylim([0, 5])
legend(fontsize=10, ncol=3)
grid()
xlabel("Position (μm)")
ylabel("Difference (V)")
NaCsPlot.maybe_save("$(prefix)_diff")

NaCsPlot.maybe_show()
