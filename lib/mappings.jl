#!/usr/bin/julia

module Mappings

export dsub100_to_dsub25, dsub100_to_electrode,
    dsub25_to_electrode, electrode_to_dsub25,
    ao_to_electrode, electrode_to_ao

const dsub100_to_dsub25 = Dict("1-1"=>"4-2-4",
                               "1-2"=>"1-2-4",
                               "1-3"=>"3-2-5",
                               "1-4"=>"1-2-5",
                               "1-5"=>"3-1-6",
                               "1-6"=>"4-1-6",
                               "1-7"=>"3-2-6",
                               "1-8"=>"4-2-2",
                               "1-9"=>"1-1-7",
                               "1-10"=>"3-2-3",
                               "1-11"=>"GND",
                               "1-12"=>"2-2-6",
                               "1-13"=>"3-2-7",
                               "1-14"=>"4-2-1",
                               "1-15"=>"2-2-7",
                               "1-16"=>"4-1-8",
                               "1-17"=>"3-1-11",
                               "1-18"=>"4-1-10",
                               "1-19"=>"3-2-11",
                               "1-20"=>"3-1-9",
                               "1-21"=>"1-2-9",
                               "1-22"=>"4-1-9",
                               "1-23"=>"1-2-10",
                               "1-24"=>"GND",
                               "1-25"=>"1-1-11",
                               "2-1"=>"4-1-4",
                               "2-2"=>"1-1-4",
                               "2-3"=>"1-1-5",
                               "2-4"=>"4-1-5",
                               "2-5"=>"3-1-7",
                               "2-6"=>"4-1-7",
                               "2-7"=>"3-2-2",
                               "2-8"=>"4-2-5",
                               "2-9"=>"3-1-3",
                               "2-10"=>"3-1-4",
                               "2-11"=>"3-1-8",
                               "2-12"=>"3-1-10",
                               "2-13"=>"4-2-7",
                               "2-14"=>"3-2-10",
                               "2-15"=>"1-2-8",
                               "2-16"=>"4-2-9",
                               "2-17"=>"3-2-8",
                               "2-18"=>"4-2-8",
                               "2-19"=>"3-1-12",
                               "2-20"=>"1-1-9",
                               "2-21"=>"3-1-13",
                               "2-22"=>"4-1-13",
                               "2-23"=>"1-1-10",
                               "2-24"=>"4-2-10",
                               "3-1"=>"1-1-3",
                               "3-2"=>"3-1-5",
                               "3-3"=>"1-2-3",
                               "3-4"=>"2-1-4",
                               "3-5"=>"2-2-5",
                               "3-6"=>"1-1-6",
                               "3-7"=>"GND",
                               "3-8"=>"1-2-6",
                               "3-9"=>"2-1-6",
                               "3-10"=>"2-1-2",
                               "3-11"=>"2-2-1",
                               "3-12"=>"4-2-11",
                               "3-13"=>"2-1-7",
                               "3-14"=>"4-2-12",
                               "3-15"=>"1-1-8",
                               "3-16"=>"1-2-7",
                               "3-17"=>"3-1-2",
                               "3-18"=>"4-1-12",
                               "3-19"=>"4-2-6",
                               "3-20"=>"2-2-11",
                               "3-21"=>"2-1-9",
                               "3-22"=>"1-1-13",
                               "3-23"=>"2-1-8",
                               "3-24"=>"2-2-8",
                               "3-25"=>"3-2-9",
                               "4-1"=>"1-2-2",
                               "4-2"=>"3-2-4",
                               "4-3"=>"2-1-5",
                               "4-4"=>"2-2-4",
                               "4-5"=>"1-1-2",
                               "4-6"=>"2-1-3",
                               "4-7"=>"4-2-3",
                               "4-8"=>"2-2-3",
                               "4-9"=>"1-2-1",
                               "4-10"=>"2-2-2",
                               "4-11"=>"4-1-2",
                               "4-12"=>"GND",
                               "4-13"=>"4-1-3",
                               "4-14"=>"2-2-10",
                               "4-15"=>"2-1-13",
                               "4-16"=>"2-2-12",
                               "4-17"=>"2-1-12",
                               "4-18"=>"3-2-1",
                               "4-19"=>"2-1-11",
                               "4-20"=>"4-1-11",
                               "4-21"=>"1-2-12",
                               "4-22"=>"2-1-10",
                               "4-23"=>"2-2-9",
                               "4-24"=>"1-1-12",
                               "4-25"=>"3-2-12",
                               "4-26"=>"1-2-11")

for (k, v) in dsub100_to_dsub25
    if v == "GND"
        continue
    end
    vs = parse.(Int, split(v, "-"))
    if vs[2] == 1
        dp = 14 - vs[3]
    else
        @assert(vs[2] == 2)
        dp = 26 - vs[3]
    end
    dsub100_to_dsub25[k] = "$(vs[1])-$(dp)"
end

const dsub100_to_electrode = Dict("4-1"=>"Q16",
                                  "4-2"=>"Q46",
                                  "4-3"=>"Q62",
                                  "4-4"=>"Q30",
                                  "4-5"=>"Q24",
                                  "4-6"=>"Q34",
                                  "4-7"=>"S6",
                                  "4-8"=>"Q36",
                                  "4-9"=>"Q56",
                                  "4-10"=>"Q28",
                                  "4-11"=>"S0",
                                  "4-12"=>"W1-1",
                                  "4-13"=>"O1",
                                  "4-14"=>"Q61",
                                  "4-15"=>"Q37",
                                  "4-16"=>"Q59",
                                  "4-17"=>"Q29",
                                  "4-18"=>"W2-1",
                                  "4-19"=>"Q33",
                                  "4-20"=>"S7",
                                  "4-21"=>"Q27",
                                  "4-22"=>"Q39",
                                  "4-23"=>"Q41",
                                  "4-24"=>"Q23",
                                  "4-25"=>"Q9",
                                  "4-26"=>"Q55",
                                  "3-1"=>"Q54",
                                  "3-2"=>"Q8",
                                  "3-3"=>"Q22",
                                  "3-4"=>"Q40",
                                  "3-5"=>"Q38",
                                  "3-6"=>"Q26",
                                  "3-7"=>"W1-2",
                                  "3-8"=>"Q20",
                                  "3-9"=>"Q32",
                                  "3-10"=>"Q60",
                                  "3-11"=>"Q42",
                                  "3-12"=>"S4",
                                  "3-13"=>"Q58",
                                  "3-14"=>"S2",
                                  "3-15"=>"Q57",
                                  "3-16"=>"Q21",
                                  "3-17"=>"W2-2",
                                  "3-18"=>"S1",
                                  "3-19"=>"S5",
                                  "3-20"=>"Q43",
                                  "3-21"=>"Q35",
                                  "3-22"=>"Q25",
                                  "3-23"=>"Q31",
                                  "3-24"=>"Q63",
                                  "3-25"=>"Q47",
                                  "2-1"=>"S8",
                                  "2-2"=>"Q50",
                                  "2-3"=>"Q18",
                                  "2-4"=>"L0",
                                  "2-5"=>"Q6",
                                  "2-6"=>"L6",
                                  "2-7"=>"Q0",
                                  "2-8"=>"L4",
                                  "2-9"=>"Q44",
                                  "2-10"=>"TV1",
                                  "2-11"=>"TI2",
                                  "2-12"=>"TV2",
                                  "2-13"=>"O0",
                                  "2-14"=>"Q45",
                                  "2-15"=>"Q17",
                                  "2-16"=>"L5",
                                  "2-17"=>"Q7",
                                  "2-18"=>"L1",
                                  "2-19"=>"Q5",
                                  "2-20"=>"Q53",
                                  "2-21"=>"Q49",
                                  "2-22"=>"L3",
                                  "2-23"=>"Q15",
                                  "2-24"=>"S11",
                                  "1-1"=>"L2",
                                  "1-2"=>"Q52",
                                  "1-3"=>"Q48",
                                  "1-4"=>"Q14",
                                  "1-5"=>"Q2",
                                  "1-6"=>"S10",
                                  "1-7"=>"Q10",
                                  "1-8"=>"L8",
                                  "1-9"=>"Q12",
                                  "1-10"=>"Q4",
                                  "1-11"=>"GND",
                                  "1-12"=>"Q64",
                                  "1-13"=>"TI1",
                                  "1-14"=>"S3",
                                  "1-15"=>"Q65",
                                  "1-16"=>"S9",
                                  "1-17"=>"Q1",
                                  "1-18"=>"L9",
                                  "1-19"=>"Q3",
                                  "1-20"=>"Q11",
                                  "1-21"=>"Q51",
                                  "1-22"=>"L7",
                                  "1-23"=>"Q19",
                                  "1-24"=>"GND",
                                  "1-25"=>"Q13")

const dsub25_to_electrode = Dict{String,String}()
const electrode_to_dsub25 = Dict{String,String}()

const ao_to_electrode = Dict{Int,String}()
const electrode_to_ao = Dict{String,Int}()

for (k, v) in dsub100_to_dsub25
    if v == "GND"
        continue
    end
    electrode = dsub100_to_electrode[k]
    if electrode == "GND"
        continue
    end
    vs = parse.(Int, split(v, "-"))
    ao = (vs[1] - 1) * 25 + vs[2] - 1
    dsub25_to_electrode[v] = electrode
    electrode_to_dsub25[electrode] = v
    ao_to_electrode[ao] = electrode
    electrode_to_ao[electrode] = ao
end

# X position of electrodes
# We need to know the axial (X) positions of the electrodes so that we can figure out
# which electrodes to use for generating potentials at a given location.

# From looking at the Sandia 3D model, each inner electrode is 70um wide in X direction
# (67um + 3um gap) and each outer quantum is 2x this (140um total).
# All the odd electrode are always located at the same X position as
# the electrode numbered one less than it.

# In unit of 70um and showing only the even electrodes, the order/positions
# of the electrodes are,

# Outer: |                47(O0)              |22(Q44-64)|          14(O0)          |
# Inner: |2(gap)|10(GND)|5(L0-8)|30(S10-0 x 5)|22(Q0-42) |8(S0-10,0-3)|6(GND)|2(gap)|

# where the number outside the parenthesis is the width in unit of 70um
# and the parenthesis marks the corresponding (even) electrode.
# The origin is located in the middle of the quantum region (11 from left and right).
# This distribution is cross-checked with the potential data
# by setting two pairs of electrode to opposite values and measuring the position
# of the zero crossing.

struct ElectrodePosition
    name::String
    left::Float64
    right::Float64
    up::Bool
end

function distance(pos::ElectrodePosition, x)
    if x < pos.left
        return pos.left - x
    elseif x > pos.right
        return x - pos.right
    else
        return 0.0
    end
end

const outer_positions = ElectrodePosition[]
const inner_positions = ElectrodePosition[]

function __populate_positions()
    begin_gnd = 12
    nL = 5
    nS = 6
    S_rep1 = 5
    nQ = 22
    S_rep2 = 1
    end_gnd = 8

    unit_um = 70

    @assert nQ % 2 == 0
    nQ_outer = nQ ÷ 2
    left_edge = -(begin_gnd + nL + nS * S_rep1 + nQ ÷ 2)

    pos_inner = left_edge + begin_gnd
    pos_outer = left_edge

    # Loading
    for i in 1:nL
        push!(inner_positions,
              ElectrodePosition("L$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(inner_positions,
              ElectrodePosition("L$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Transition 1
    for j in 1:S_rep1
        for i in nS:-1:1
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, true))
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, false))
            pos_inner += 1
        end
    end

    # Outer 1
    push!(outer_positions,
          ElectrodePosition("O0", pos_outer * unit_um,
                            pos_inner * unit_um, true))
    push!(outer_positions,
          ElectrodePosition("O1", pos_outer * unit_um,
                            pos_inner * unit_um, false))
    pos_outer = pos_inner

    # Quantum inner
    for i in 1:nQ
        push!(inner_positions,
              ElectrodePosition("Q$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(inner_positions,
              ElectrodePosition("Q$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Quantum outer
    for i in 1:nQ_outer
        i += nQ
        push!(outer_positions,
              ElectrodePosition("Q$(i * 2 - 2)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um, true))
        push!(outer_positions,
              ElectrodePosition("Q$(i * 2 - 1)", pos_outer * unit_um,
                                (pos_outer + 2) * unit_um, false))
        pos_outer += 2
    end
    @assert pos_inner == pos_outer

    # Transition 2
    for j in 1:S_rep2
        for i in 1:nS
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, true))
            push!(inner_positions,
                  ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                    (pos_inner + 1) * unit_um, false))
            pos_inner += 1
        end
    end

    # S0-S3 appeared again at the end (shouldn't really matter......)
    for i in 1:4
        push!(inner_positions,
              ElectrodePosition("S$(i * 2 - 2)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, true))
        push!(inner_positions,
              ElectrodePosition("S$(i * 2 - 1)", pos_inner * unit_um,
                                (pos_inner + 1) * unit_um, false))
        pos_inner += 1
    end

    # Outer 2
    push!(outer_positions,
          ElectrodePosition("O0", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um, true))
    push!(outer_positions,
          ElectrodePosition("O1", pos_outer * unit_um,
                            (pos_inner + end_gnd) * unit_um, false))
end

__populate_positions()

mutable struct ElectrodeSearchState
    inner_idx1::Int
    inner_idx2::Int
    outer_idx1::Int
    outer_idx2::Int

    const pos::Float64
    const inner_candidates::Vector{ElectrodePosition}
    const outer_candidates::Vector{ElectrodePosition}
    function ElectrodeSearchState(pos)
        inner_idx2 = searchsortedfirst(inner_positions, pos, lt=(x, y)->x.right < y)
        outer_idx2 = searchsortedfirst(outer_positions, pos, lt=(x, y)->x.right < y)

        return new(inner_idx2 - 1, inner_idx2, outer_idx2 - 1, outer_idx2,
                   pos, ElectrodePosition[], ElectrodePosition[])
    end
end

function _find_next_distance!(state::ElectrodeSearchState)
    # First find the closest distance
    dist = Inf
    if state.inner_idx2 <= length(inner_positions)
        dist = min(dist, distance(inner_positions[state.inner_idx2], state.pos))
    end
    if state.inner_idx1 > 0
        dist = min(dist, distance(inner_positions[state.inner_idx1], state.pos))
    end
    if state.outer_idx2 <= length(outer_positions)
        dist = min(dist, distance(outer_positions[state.outer_idx2], state.pos))
    end
    if state.outer_idx1 > 0
        dist = min(dist, distance(outer_positions[state.outer_idx1], state.pos))
    end
    if !isfinite(dist)
        return dist
    end
    empty!(state.inner_candidates)
    empty!(state.outer_candidates)
    while state.inner_idx2 <= length(inner_positions)
        epos = inner_positions[state.inner_idx2]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.inner_candidates, epos)
        state.inner_idx2 += 1
    end
    while state.inner_idx1 > 0
        epos = inner_positions[state.inner_idx1]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.inner_candidates, epos)
        state.inner_idx1 -= 1
    end
    while state.outer_idx2 <= length(outer_positions)
        epos = outer_positions[state.outer_idx2]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.outer_candidates, epos)
        state.outer_idx2 += 1
    end
    while state.outer_idx1 > 0
        epos = outer_positions[state.outer_idx1]
        if distance(epos, state.pos) > dist
            break
        end
        push!(state.outer_candidates, epos)
        state.outer_idx1 -= 1
    end
    @assert !isempty(state.inner_candidates) || !isempty(state.outer_candidates)
    return dist
end

"""
Find at least `min_num` electrodes that are the closest in axial (X) position
to `pos` (in um). All electrodes within `min_dist` will also be included.

`electrode_index` is a map from electrode name to a unique ID.
ID 1 will be ignored (assumed to be ground)
electrodes with the same ID are assumed to be shorted together
and therefore will be treated as the same one.
"""
function find_electrodes(electrode_index, pos; min_num=0, min_dist=0)
    res = Set{Int}()

    search_state = ElectrodeSearchState(pos)
    dist_satisfied = false
    num_satisfied = false

    while true
        num_satisfied = min_num <= length(res)
        if num_satisfied && dist_satisfied
            return res
        end

        dist = _find_next_distance!(search_state)
        if !isfinite(dist)
            if num_satisfied
                return res
            end
            error("Unable to find enough terms")
        end
        dist_satisfied = dist >= min_dist

        for p in search_state.inner_candidates
            id = electrode_index[p.name]
            # Ground
            if id == 1
                continue
            end
            push!(res, id)
        end

        for p in search_state.outer_candidates
            id = electrode_index[p.name]
            # Ground
            if id == 1
                continue
            end
            push!(res, id)
        end
    end
end

end
