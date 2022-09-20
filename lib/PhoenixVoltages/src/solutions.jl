#!/usr/bin/julia

"""
Tools to find voltage solutions (in terms of voltages on the electrodes)
"""
module Solutions

import ..Mappings
import ..Fitting
using ..Potentials
import ..gradient
using ..Outputs
using ..Optimizers
using NLsolve
using LinearAlgebra
using DelimitedFiles
using MAT

function find_flat_point(data::A; init=ntuple(i->(size(data, i) + 1) / 2, Val(N))) where (A<:AbstractArray{T,N} where T) where N
    fitter = Fitting.PolyFitter(ntuple(i->3, Val(N))...)
    cache = Fitting.PolyFitCache(fitter, data)
    function model!(g, x)
        xt = ntuple(i->x[i], Val(N))
        for i in 1:N
            g[i] = gradient(cache, i, xt...)
        end
    end
    res = nlsolve(model!, collect(init))
    return ntuple(i->res.zero[i], Val(N))
end

function find_all_flat_points(all_data::A; init=ntuple(i->(size(all_data, i) + 1) / 2, Val(N - 1))) where (A<:AbstractArray{T,N} where T) where N

    npoints = size(all_data, 3)
    all_res = Matrix{Float64}(undef, npoints, N - 1)

    for i in 1:npoints
        init = find_flat_point(@view(all_data[:, :, i]), init=init)
        all_res[i, :] .= init
    end
    return all_res
end

# EURIQA unit:
# Unit such that electric potential that creates 1MHz trapping frequency
# for Yb171 has the form X^2/2,
# and the electric potential between two ions is 1/r.

const N_A = 6.02214076e23
const m_Yb171 = 170.9363315e-3 / N_A # kg
const q_e = 1.60217663e-19 # C
const ε_0 = 8.8541878128e-12

let A = m_Yb171 * (2π * 1e6)^2 / q_e,
    B = q_e / (4π * ε_0)

    global const V_unit = cbrt(A * B^2) # V
    global const l_unit = cbrt(B / A) # m
end
const l_unit_um = l_unit * 1e6 # um
const V_unit_uV = V_unit * 1e6 # uV

# Terms we care about
# x, y, z, xy, yz, xz, (z^2 - y^2) / 2, (x^2 - (y^2 + z^2) / 2) / 2, x^3 / 3!, x^4 / 4!
# Since we care about the symmetry of the x^2 and z^2 term,
# we actually do need to scale the x, y and z correctly.
# stride should be in um, voltage should be in V
function get_compensate_terms1(res::Fitting.PolyFitResult{3}, stride)
    # axis order of fitting result is z, y, x
    # axis order of stride is x, y, z
    raw_x = res[0, 0, 1]
    raw_y = res[0, 1, 0]
    raw_z = res[1, 0, 0]

    raw_xy = res[0, 1, 1]
    raw_yz = res[1, 1, 0]
    raw_zx = res[1, 0, 1]

    raw_x2 = res[0, 0, 2]
    raw_y2 = res[0, 2, 0]
    raw_z2 = res[2, 0, 0]

    raw_x3 = res[0, 0, 3]
    raw_x4 = res[0, 0, 4]

    scaled_x = raw_x / stride[1]
    scaled_y = raw_y / stride[2]
    scaled_z = raw_z / stride[3]

    # We need to divide the xy/yz/zx terms by 2 relative to the x2, y2, z2 terms.
    # This makes sure that, e.g., the z^2 term is a direct rotation of the
    # xy/yz/zx terms.
    scaled_xy = raw_xy / stride[1] / stride[2]
    scaled_yz = raw_yz / stride[2] / stride[3]
    scaled_zx = raw_zx / stride[3] / stride[1]

    scaled_x2 = raw_x2 / stride[1]^2 * 2
    scaled_y2 = raw_y2 / stride[2]^2 * 2
    scaled_z2 = raw_z2 / stride[3]^2 * 2

    scaled_x3 = raw_x3 / stride[1]^3 * 6
    scaled_x4 = raw_x4 / stride[1]^4 * 24

    # The two legal quadratic terms are `x^2 - (y^2 + z^2) / 2` and `z^2 - y^2`
    # which are also orthogonal to each other.
    # The orthogonal illegal term is `x^2 + y^2 + z^2`.
    # Here we just need to find the transfermation to go from the taylor expansion
    # basis to the new basis.
    # Since the three terms are orthogonal, we can just compute the dot product
    # with these three terms and apply the correct normalization coefficient.
    xx = (2 * scaled_x2 - scaled_y2 - scaled_z2) / 3
    zz = (scaled_z2 - scaled_y2) / 2

    # Current units are V/um^n
    # Expected units
    # DX/DY/DZ: V/m
    # XY, YZ, ZX, ZZ, XX: 525 uV / (2.74 um)^2
    # X3: 525 uV / (2.74 um)^3
    # X4: 525 uV / (2.74 um)^4
    scale_1 = 1e6
    scale_2 = (l_unit_um^2 / V_unit)
    scale_3 = (l_unit_um^3 / V_unit)
    scale_4 = (l_unit_um^4 / V_unit)
    return (dx=scaled_x * scale_1, dy=scaled_y * scale_1, dz=scaled_z * scale_1,
            xy=scaled_xy * scale_2, yz=scaled_yz * scale_2, zx=scaled_zx * scale_2,
            z2=zz * scale_2, x2=xx * scale_2, x3=scaled_x3 * scale_3,
            x4=scaled_x4 * scale_4)
end

function compensate_fitter1(solution::Potential; sizes=(5, 5, 129))
    fitter = Fitting.PolyFitter(2, 2, 4, sizes=sizes)
    return Potentials.FitCache(fitter, solution)
end

function get_compensate_coeff1(cache::Potentials.FitCache, pos::NTuple{3})
    # pos is in xyz index

    x_coord = x_index_to_axis(cache.solution, pos[1]) .* 1000
    ele_select = Mappings.find_electrodes(cache.solution.electrode_index,
                                          x_coord, min_num=20, min_dist=350)
    ele_select = sort!(collect(ele_select))
    fits = [get(cache, e, (pos[3], pos[2], pos[1])) for e in ele_select]

    # Change stride to um in unit
    stride_um = cache.solution.stride .* 1000
    nfits = length(fits)
    coefficient = Matrix{Float64}(undef, 10, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(get_compensate_terms1(fits[i], stride_um))
    end
    return ele_select, coefficient
end

function solve_compensate1(cache::Potentials.FitCache, pos::NTuple{3})
    ele_select, coefficient = get_compensate_coeff1(cache, pos)
    # X = coefficient \ Matrix(I, 10, 10)
    X = Optimizers.optimize_minmax(coefficient, Matrix(I, 10, 10))
    @assert size(X, 2) == 10
    return ele_select, (dx=X[:, 1], dy=X[:, 2], dz=X[:, 3],
                        xy=X[:, 4], yz=X[:, 5], zx=X[:, 6],
                        z2=X[:, 7], x2=X[:, 8], x3=X[:, 9], x4=X[:, 10])
end

# Terms we care about during transport
# x, y, z, xy, yz, (z^2 - y^2) / 2, (x^2 - (y^2 + z^2) / 2) / 2, x^3 / 3!, x^4 / 4!
# zx is missing here since it seems to require a fairly high voltage to compensate
# Since we care about the symmetry of the x^2 and z^2 term,
# we actually do need to scale the x, y and z correctly.
# stride should be in um, voltage should be in V
function get_compensate_terms1_nozx(res::Fitting.PolyFitResult{3}, stride)
    # axis order of fitting result is z, y, x
    # axis order of stride is x, y, z
    raw_x = res[0, 0, 1]
    raw_y = res[0, 1, 0]
    raw_z = res[1, 0, 0]

    raw_xy = res[0, 1, 1]
    raw_yz = res[1, 1, 0]

    raw_x2 = res[0, 0, 2]
    raw_y2 = res[0, 2, 0]
    raw_z2 = res[2, 0, 0]

    raw_x3 = res[0, 0, 3]
    raw_x4 = res[0, 0, 4]

    scaled_x = raw_x / stride[1]
    scaled_y = raw_y / stride[2]
    scaled_z = raw_z / stride[3]

    # We need to divide the xy/yz/zx terms by 2 relative to the x2, y2, z2 terms.
    # This makes sure that, e.g., the z^2 term is a direct rotation of the
    # xy/yz/zx terms.
    scaled_xy = raw_xy / stride[1] / stride[2]
    scaled_yz = raw_yz / stride[2] / stride[3]

    scaled_x2 = raw_x2 / stride[1]^2 * 2
    scaled_y2 = raw_y2 / stride[2]^2 * 2
    scaled_z2 = raw_z2 / stride[3]^2 * 2

    scaled_x3 = raw_x3 / stride[1]^3 * 6
    scaled_x4 = raw_x4 / stride[1]^4 * 24

    # The two legal quadratic terms are `x^2 - (y^2 + z^2) / 2` and `z^2 - y^2`
    # which are also orthogonal to each other.
    # The orthogonal illegal term is `x^2 + y^2 + z^2`.
    # Here we just need to find the transfermation to go from the taylor expansion
    # basis to the new basis.
    # Since the three terms are orthogonal, we can just compute the dot product
    # with these three terms and apply the correct normalization coefficient.
    xx = (2 * scaled_x2 - scaled_y2 - scaled_z2) / 3
    zz = (scaled_z2 - scaled_y2) / 2

    # Current units are V/um^n
    # Expected units
    # DX/DY/DZ: V/m
    # XY, YZ, ZZ, XX: 525 uV / (2.74 um)^2
    # X3: 525 uV / (2.74 um)^3
    # X4: 525 uV / (2.74 um)^4
    scale_1 = 1e6
    scale_2 = (l_unit_um^2 / V_unit)
    scale_3 = (l_unit_um^3 / V_unit)
    scale_4 = (l_unit_um^4 / V_unit)
    return (dx=scaled_x * scale_1, dy=scaled_y * scale_1, dz=scaled_z * scale_1,
            xy=scaled_xy * scale_2, yz=scaled_yz * scale_2,
            z2=zz * scale_2, x2=xx * scale_2, x3=scaled_x3 * scale_3,
            x4=scaled_x4 * scale_4)
end

function get_compensate_coeff1_nozx(cache::Potentials.FitCache, pos::NTuple{3})
    # pos is in xyz index

    x_coord = x_index_to_axis(cache.solution, pos[1]) .* 1000
    ele_select = Mappings.find_electrodes(cache.solution.electrode_index,
                                          x_coord, min_num=20, min_dist=350)
    ele_select = sort!(collect(ele_select))
    fits = [get(cache, e, (pos[3], pos[2], pos[1])) for e in ele_select]

    # Change stride to um in unit
    stride_um = cache.solution.stride .* 1000
    nfits = length(fits)
    coefficient = Matrix{Float64}(undef, 9, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(get_compensate_terms1_nozx(fits[i], stride_um))
    end
    return ele_select, coefficient
end

function solve_compensate1_nozx(cache::Potentials.FitCache, pos::NTuple{3})
    ele_select, coefficient = get_compensate_coeff1_nozx(cache, pos)
    # X = coefficient \ Matrix(I, 9, 9)
    X = Optimizers.optimize_minmax(coefficient, Matrix(I, 9, 9))
    @assert size(X, 2) == 9
    return ele_select, (dx=X[:, 1], dy=X[:, 2], dz=X[:, 3],
                        xy=X[:, 4], yz=X[:, 5],
                        z2=X[:, 6], x2=X[:, 7], x3=X[:, 8], x4=X[:, 9])
end

function solve_transfer1(cache::Potentials.FitCache, pos::NTuple{3})
    ele_select, coefficient = get_compensate_coeff1_nozx(cache, pos)
    y = zeros(9)
    y[7] = 1 # X2
    y[5] = 1 # YZ
    # return ele_select, coefficient \ y
    return ele_select, Optimizers.optimize_minmax(coefficient, y)
end

struct CenterTracker
    zy_index::Matrix{Float64}
end

function _try_suffix(name)
    isfile(name) && return name
    isfile(name * ".mat") && return name * ".mat"
    return
end

function _try_path(name)
    if !isabspath(name)
        p = _try_suffix(joinpath(@__DIR__, "../data", name))
        p !== nothing && return p
    end
    return _try_suffix(name)
end

function CenterTracker(name::AbstractString="rf_center")
    name = _try_path(name)
    return matopen(name) do mat
        return CenterTracker(read(mat, "zy_index"))
    end
end

function Base.get(tracker::CenterTracker, xidx)
    # return (y, z)
    nx = size(tracker.zy_index, 1)
    lb_idx = min(max(floor(Int, xidx), 1), nx)
    ub_idx = min(max(ceil(Int, xidx), 1), nx)
    y_lb = tracker.zy_index[lb_idx, 2]
    z_lb = tracker.zy_index[lb_idx, 1]
    if lb_idx == ub_idx
        return y_lb, z_lb
    end
    @assert ub_idx == lb_idx + 1
    y_ub = tracker.zy_index[ub_idx, 2]
    z_ub = tracker.zy_index[ub_idx, 1]
    c_ub = xidx - lb_idx
    c_lb = ub_idx - xidx
    return y_lb * c_lb + y_ub * c_ub, z_lb * c_lb + z_ub * c_ub
end

function load_short_map(fname)
    m = readdlm(fname, ',', String)
    res = Dict{String,String}()
    for i in 1:size(m, 1)
        res[m[i, 1]] = m[i, 2]
    end
    return res
end

function fill_data_line!(values, solution::Potential, mapfile::MapFile,
                         electrodes, term)
    term_map = Dict(zip(electrodes, term))
    nelectrodes = length(mapfile.names)
    @assert length(values) == nelectrodes
    for i in 1:nelectrodes
        id = get(solution.electrode_index, mapfile.names[i], -1)
        values[i] = get(term_map, id, 0.0)
    end
    return values
end

function get_data_line(solution::Potential, mapfile::MapFile,
                       electrodes, term)
    nelectrodes = length(mapfile.names)
    return fill_data_line!(Vector{Float64}(undef, nelectrodes),
                           solution, mapfile, electrodes, term)
end

function compensation_to_file(solution::Potential, mapfile::MapFile,
                              electrodes, terms)
    term_names = String[]
    term_values = Vector{Float64}[]
    nelectrodes = length(mapfile.names)
    # For DX, DY, DZ the EURIQA frontend expects a different unit in the config file
    # compared to the UI...
    # DX
    push!(term_names, "DX")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.dx) .* 1000)
    # DY
    push!(term_names, "DY")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.dy) .* 1000)
    # DZ
    push!(term_names, "DZ")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.dz) .* 1000)
    # QZY
    push!(term_names, "QZY")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.yz))
    # QZZ
    push!(term_names, "QZZ")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.z2))
    # QXZ
    push!(term_names, "QXZ")
    if !hasproperty(terms, :zx)
        push!(term_values, zeros(nelectrodes))
    else
        push!(term_values, get_data_line(solution, mapfile, electrodes, terms.zx))
    end
    # X1
    # DX is in V/m, X1 is in 525 uV / 2.74 um
    push!(term_names, "X1")
    push!(term_values, get_data_line(solution, mapfile, electrodes,
                                     terms.dx .* (V_unit_uV / l_unit_um)))
    # X2
    push!(term_names, "X2")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.x2))
    # X3
    push!(term_names, "X3")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.x3))
    # X4
    push!(term_names, "X4")
    push!(term_values, get_data_line(solution, mapfile, electrodes, terms.x4))
    # JunctionCenter
    push!(term_names, "JunctionCenter")
    push!(term_values, zeros(nelectrodes))
    # JunctionTransition
    push!(term_names, "JunctionTransition")
    push!(term_values, zeros(nelectrodes))
    # LoadInners
    push!(term_names, "LoadInners")
    push!(term_values, zeros(nelectrodes))
    # LoadOuters
    push!(term_names, "LoadOuters")
    push!(term_values, zeros(nelectrodes))
    # Eject
    push!(term_names, "Eject")
    push!(term_values, zeros(nelectrodes))

    return CompensationFile(mapfile, term_names, term_values)
end

# Terms we care about on pheonix
# x, y, z, xy, yz, xz, (z^2 - y^2) / 2, (x^2 - (y^2 + z^2) / 2) / 2,
# x^3 / 3!, x^4 / 4!, x^2z / 2
# Compared to HOA, we are able to compensate for x^2z due to the outer electrodes.
# stride should be in um, voltage should be in V
function get_compensate_terms2(res::Fitting.PolyFitResult{3}, stride)
    # axis order of fitting result is z, y, x
    # axis order of stride is x, y, z
    raw_x = res[0, 0, 1]
    raw_y = res[0, 1, 0]
    raw_z = res[1, 0, 0]

    raw_xy = res[0, 1, 1]
    raw_yz = res[1, 1, 0]
    raw_zx = res[1, 0, 1]

    raw_x2 = res[0, 0, 2]
    raw_y2 = res[0, 2, 0]
    raw_z2 = res[2, 0, 0]

    raw_x3 = res[0, 0, 3]
    raw_x4 = res[0, 0, 4]

    raw_x2z = res[1, 0, 2]

    scaled_x = raw_x / stride[1]
    scaled_y = raw_y / stride[2]
    scaled_z = raw_z / stride[3]

    # We need to divide the xy/yz/zx terms by 2 relative to the x2, y2, z2 terms.
    # This makes sure that, e.g., the z^2 term is a direct rotation of the
    # xy/yz/zx terms.
    scaled_xy = raw_xy / stride[1] / stride[2]
    scaled_yz = raw_yz / stride[2] / stride[3]
    scaled_zx = raw_zx / stride[3] / stride[1]

    scaled_x2 = raw_x2 / stride[1]^2 * 2
    scaled_y2 = raw_y2 / stride[2]^2 * 2
    scaled_z2 = raw_z2 / stride[3]^2 * 2

    scaled_x3 = raw_x3 / stride[1]^3 * 6
    scaled_x4 = raw_x4 / stride[1]^4 * 24

    scaled_x2z = raw_x2z / stride[1]^2 / stride[3] * 2

    # The two legal quadratic terms are `x^2 - (y^2 + z^2) / 2` and `z^2 - y^2`
    # which are also orthogonal to each other.
    # The orthogonal illegal term is `x^2 + y^2 + z^2`.
    # Here we just need to find the transfermation to go from the taylor expansion
    # basis to the new basis.
    # Since the three terms are orthogonal, we can just compute the dot product
    # with these three terms and apply the correct normalization coefficient.
    xx = (2 * scaled_x2 - scaled_y2 - scaled_z2) / 3
    zz = (scaled_z2 - scaled_y2) / 2

    # Current units are V/um^n
    # Expected units
    # DX/DY/DZ: V/m
    # XY, YZ, ZX, ZZ, XX: 525 uV / (2.74 um)^2
    # X3, X2Z: 525 uV / (2.74 um)^3
    # X4: 525 uV / (2.74 um)^4
    scale_1 = 1e6
    scale_2 = (l_unit_um^2 / V_unit)
    scale_3 = (l_unit_um^3 / V_unit)
    scale_4 = (l_unit_um^4 / V_unit)
    return (dx=scaled_x * scale_1, dy=scaled_y * scale_1, dz=scaled_z * scale_1,
            xy=scaled_xy * scale_2, yz=scaled_yz * scale_2, zx=scaled_zx * scale_2,
            z2=zz * scale_2, x2=xx * scale_2, x3=scaled_x3 * scale_3,
            x4=scaled_x4 * scale_4, x2z=scaled_x2z * scale_3)
end

function get_compensate_coeff2(cache::Potentials.FitCache, pos::NTuple{3})
    # pos is in xyz index

    x_coord = x_index_to_axis(cache.solution, pos[1]) .* 1000
    ele_select = Mappings.find_electrodes(cache.solution.electrode_index,
                                          x_coord, min_num=20, min_dist=350)
    ele_select = sort!(collect(ele_select))
    fits = [get(cache, e, (pos[3], pos[2], pos[1])) for e in ele_select]

    # Change stride to um in unit
    stride_um = cache.solution.stride .* 1000
    nfits = length(fits)
    coefficient = Matrix{Float64}(undef, 11, nfits)
    for i in 1:nfits
        coefficient[:, i] .= Tuple(get_compensate_terms2(fits[i], stride_um))
    end
    return ele_select, coefficient
end

function solve_compensate2(cache::Potentials.FitCache, pos::NTuple{3})
    ele_select, coefficient = get_compensate_coeff2(cache, pos)
    # X = coefficient \ Matrix(I, 10, 10)
    X = Optimizers.optimize_minmax(coefficient, Matrix(I, 11, 11))
    @assert size(X, 2) == 11
    return ele_select, (dx=X[:, 1], dy=X[:, 2], dz=X[:, 3],
                        xy=X[:, 4], yz=X[:, 5], zx=X[:, 6],
                        z2=X[:, 7], x2=X[:, 8], x3=X[:, 9], x4=X[:, 10],
                        x2z=X[:, 11])
end


end
