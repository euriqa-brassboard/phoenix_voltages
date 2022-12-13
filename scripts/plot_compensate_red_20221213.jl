#!/usr/bin/julia

using MAT
using PyPlot
using NaCsPlot
using Statistics

const prefixes = ["compensate_red_20221213"=>("Origin", 10)]
const prefix = joinpath(@__DIR__, "../imgs/compensate_red_20221213")

function load_solution(name)
    solution, termname = matopen(joinpath(@__DIR__, "../data", name)) do mat
        return read(mat, "transfer_solutions"), read(mat, "termname")
    end
    xpos_um = [sol["xpos_um"] for sol in solution]
    maxv = [maximum(abs.(sol["voltages"])) for sol in solution]
    return (termname=termname, xpos_um=xpos_um, maxv=maxv)
end

function load_solutions(prefix, num)
    it = (load_solution(joinpath(prefix, "$(i).mat")) for i in 1:num)
    return Dict(s.termname=>s for s in it)
end

const solutions = [(name, load_solutions(prefix, num))
                   for (prefix, (name, num)) in prefixes]

function term_name_to_latex(name)
    return '$' * replace(name, r"[0-9]*"=>s"^{\0}") * '$'
end

function make_plot(solutions, term)
    max_median = 0
    min_min = Inf
    for (i, (name, sols)) in enumerate(solutions)
        sol = get(sols, term, nothing)
        sol === nothing && continue
        plot(sol.xpos_um, sol.maxv, label=name, color="C$(i - 1)")
        max_median = max(max_median, median(sol.maxv[-750 .< sol.xpos_um .< 750]))
        min_min = min(min_min, minimum(sol.maxv))
    end
    ymax = max_median * 2.2
    ylim([max(0, min_min - (ymax - min_min) * 0.02), ymax])
    # Position of the ion (in the middle of Q16/17)
    axvline(-175, color="red", ls="--")
    # yscale("log")
    legend(fontsize=10, ncol=3)
    grid()
    grid(which="minor", alpha=0.2)
    gca().xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    gca().yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    title("$(term_name_to_latex(term))")
    xlabel("Position (μm)")
    ylabel("Max Voltage (V)")
end

figure(figsize=[6.4 * 2, 4.8 * 2])
subplot(2, 2, 1)
make_plot(solutions, "dx")
subplot(2, 2, 2)
make_plot(solutions, "dy")
subplot(2, 2, 3)
make_plot(solutions, "dz")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_1")

figure(figsize=[6.4 * 2, 4.8 * 2])
subplot(2, 2, 1)
make_plot(solutions, "xy")
subplot(2, 2, 2)
make_plot(solutions, "yz")
subplot(2, 2, 3)
make_plot(solutions, "zx")
subplot(2, 2, 4)
make_plot(solutions, "z2")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_2")

figure(figsize=[6.4 * 2, 4.8 * 2])
subplot(2, 2, 1)
make_plot(solutions, "x2")
subplot(2, 2, 2)
make_plot(solutions, "x3")
subplot(2, 2, 3)
make_plot(solutions, "x4")
# subplot(2, 2, 4)
# make_plot(solutions, "x2z")
tight_layout()
NaCsPlot.maybe_save("$(prefix)_3")

NaCsPlot.maybe_show()
