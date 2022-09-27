#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages
import PhoenixVoltages.Solutions
using PhoenixVoltages.Potentials
using PhoenixVoltages.Fitting
using PhoenixVoltages.Reconstruct
using MAT
using PyPlot
using NaCsPlot

const prefixes = [
    "compensate_20220920"=>("origin", 10),
    "compensate_x2z_20220920"=>("x2z", 11),
    "compensate_nozx_20220920"=>("nozx", 9)
]
const potential_file = ARGS[1]
const potential = Potentials.import_pillbox_64(potential_file)
const fits_cache = Solutions.compensate_fitter1(potential, sizes=(5, 5, 15))
const centers = Solutions.CenterTracker()

function get_all_fits(xpos_um, voltages, electrodes, electrode_names)
    electrodes_voltages = Dict{String,Float64}()
    for (eidx, v) in zip(electrodes, voltages)
        for name in electrode_names[eidx]
            electrodes_voltages[name] = v
        end
    end
    xidx = Potentials.x_axis_to_index(potential, xpos_um ./ 1000)
    cf = Reconstruct.get_center_fit(fits_cache, electrodes_voltages,
                                    xidx, centers=centers)
    return Solutions.get_compensate_terms2(cf, potential.stride .* 1000)
end

const term_names = [:dx, :dy, :dz, :xy, :yz, :zx, :z2, :x2, :x3, :x4, :x2z]
const term_scale = Dict{String,Float64}("dx"=>500, "dy"=>1000, "dz"=>500,
                                        "xy"=>0.25, "yz"=>0.5,
                                        "zx"=>0.125, "z2"=>0.125,
                                        "x2"=>0.1, "x3"=>0.01, "x4"=>0.001,
                                        "x2z"=>0.003)

function load_solution(name)
    data = matread(joinpath(@__DIR__, "../data", name))
    termname = data["termname"]
    solution = data["transfer_solutions"]
    xpos_um = [sol["xpos_um"] for sol in solution if -900 <= sol["xpos_um"] <= 900]
    electrode_names = data["electrode_names"]
    fits = [get_all_fits(sol["xpos_um"], sol["voltages"],
                         sol["electrodes"], electrode_names) for sol in solution
                             if -900 <= sol["xpos_um"] <= 900]
    return (termname=termname, xpos_um=xpos_um,
            fits=Dict{Symbol,Vector{Float64}}(term_name=>[getfield(fit, term_name)
                                                     for fit in fits]
                                              for term_name in term_names))
end

function load_solutions(prefix, num)
    local xpos_um
    terms = Dict{Symbol,Dict{String,Vector{Float64}}}(
        term_name=>Dict{String,Vector{Float64}}() for term_name in term_names)
    for i in 1:num
        @time sol = load_solution(joinpath(prefix, "$(i).mat"))
        if @isdefined(xpos_um)
            @assert xpos_um == sol.xpos_um
        else
            xpos_um = sol.xpos_um
        end
        for term_name in term_names
            terms[term_name][sol.termname] = sol.fits[term_name]
        end
    end
    return (xpos_um=xpos_um, terms=terms)
end

const solutions = [(name, load_solutions(prefix, num))
                   for (prefix, (name, num)) in prefixes]

for (name, solution) in solutions
    figure(figsize=[6.4 * 3, 4.8 * 4])
    for (i, term_name) in enumerate(term_names)
        subplot(4, 3, i)
        title(term_name)
        terms = solution.terms[term_name]
        base_scale = 1 / term_scale[String(term_name)]
        for term_name2 in term_names
            term_name2 = String(term_name2)
            if term_name2 in keys(terms)
                s = base_scale * term_scale[term_name2]
                plot(solution.xpos_um, terms[term_name2] .* s,
                     label=term_name2)
            end
        end
        grid()
        legend(ncol=4, fontsize=8)
    end
    tight_layout()
end

NaCsPlot.maybe_show()
