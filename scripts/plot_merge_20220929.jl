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

const solution_1 = load_solution("merge_20220929.mat")
const solution_2 = load_solution("merge_20221007_local.mat")
const solution_3 = load_solution("merge_data_20221004.mat")

const prefix = joinpath(@__DIR__, "../imgs/merge_20221001")

figure()
plot(solution_1.xpos_um, solution_1.maxv, label="0929")
plot(solution_2.xpos_um, solution_2.maxv, label="1007")
plot(solution_3.xpos_um, solution_3.maxv, label="1004")
ylim([0.5, 9])
legend(fontsize=10, ncol=3)
grid()
xlabel("Position (μm)")
ylabel("Max (V)")
NaCsPlot.maybe_save("$(prefix)_max")

figure()
plot(solution_1.xpos_um[1:end - 1], solution_1.maxdiff, label="0929")
plot(solution_2.xpos_um[1:end - 1], solution_2.maxdiff, label="1007")
plot(solution_3.xpos_um[1:end - 1], solution_3.maxdiff, label="1004")
ylim([0, 2.5])
legend(fontsize=10, ncol=3)
grid()
xlabel("Position (μm)")
ylabel("Difference (V)")
NaCsPlot.maybe_save("$(prefix)_diff")

NaCsPlot.maybe_show()
