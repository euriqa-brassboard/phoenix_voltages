#!/usr/bin/julia

module Reconstruct

import ..Fitting
import ..Potentials
import ..Solutions

function _normalize_voltage_map(potential::Potentials.Potential,
                                electrodes_voltages::AbstractDict{I} where I<:Integer)
    return Dict{Int,Float64}(electrodes_voltages)
end

function _normalize_voltage_map(potential::Potentials.Potential,
                                electrodes_voltages)
    res = Dict{Int,Float64}()
    for (ele, v) in electrodes_voltages
        if !isfinite(v)
            throw(ArgumentError("Voltage specified must be finite."))
        end
        if !isa(ele, Integer)
            ele = potential.electrode_index[ele]
        end
        if get!(res, ele, v) != v
            throw(ArgumentError("Multiple voltages specified for the same electrode."))
        end
    end
    return res
end

"""
Calculate the 3D potential given the electrode voltages.
The potential is provided in a `Potential` object whereas the voltages are given
in a map between electrode name/index and the voltage value.
"""
function get_potential(potential::Potentials.Potential, electrodes_voltages)
    electrodes_voltages = _normalize_voltage_map(potential, electrodes_voltages)
    res = zeros(size(potential.data, 1), size(potential.data, 2),
                size(potential.data, 3))
    for (ele_id, v) in electrodes_voltages
        res .+= v .* @view(potential.data[:, :, :, ele_id])
    end
    return res
end

function get_fits(fitter::Fitting.PolyFitter{3}, data::AbstractArray{T,3} where T)
    return Fitting.PolyFitCache(fitter, data)
end

function get_center_fit(fitcache::Fitting.PolyFitCache, xidx;
                        centers=Solutions.CenterTracker())
    (yidx, zidx) = get(centers, xidx)
    return get(fitcache, (zidx, yidx, xidx))
end

end
