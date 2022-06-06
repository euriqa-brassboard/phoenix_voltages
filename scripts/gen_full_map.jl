#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using DelimitedFiles

using PhoenixVoltages.Mappings

const grid_to_electrode = readdlm(joinpath(@__DIR__, "../data/grid_to_electrode.csv"), ',')
const inter_pad_to_grid = readdlm(joinpath(@__DIR__, "../data/inter_pad_to_grid.csv"), ',')

function get_d100_pin(v::String)
    vs = parse.(Int, split(v, "-"))
    @assert length(vs) == 2
    if vs[1] == 1
        return 101 - vs[2]
    elseif vs[1] == 2
        return 101 - (vs[2] + 25)
    elseif vs[1] == 3
        return 101 - (vs[2] + 25 + 24)
    else
        @assert vs[1] == 4
        return 101 - (vs[2] + 25 + 24 + 25)
    end
end

mutable struct Net
    d100::String # using the DC test convention
    d25_plug::Int
    d25_pin::Int
    dac::Int
    electrode::String
    inter_pad::Int
    grid::String
    note::String
    Net() = new("", -1, -1, -1, "", -1, "", "")
end

function print_csv(io::IO, net::Net)
    if !isempty(net.d100)
        print(io, replace(net.d100, "-"=>","), ",", get_d100_pin(net.d100), ",")
    else
        print(io, ",,,")
    end
    if net.d25_plug != -1
        print(io, net.d25_plug, ",", net.d25_pin, ",")
    else
        print(io, ",,")
    end
    if net.dac != -1
        print(io, net.dac)
    end
    print(io, ",")
    print(io, net.electrode, ",")
    if net.inter_pad != -1
        print(io, net.inter_pad)
    end
    print(io, ",")
    print(io, net.grid, ",", net.note)
end

const all_nets = Net[]

function get_net(field::Symbol, val)
    for net in all_nets
        if getfield(net, field) == val
            return net
        end
    end
    net = Net()
    setproperty!(net, field, val)
    push!(all_nets, net)
    return net
end

function get_net(field::Symbol, val, field2::Symbol, val2)
    for net in all_nets
        if getfield(net, field) == val && getfield(net, field2) == val2
            return net
        end
    end
    net = Net()
    setproperty!(net, field, val)
    setproperty!(net, field2, val2)
    push!(all_nets, net)
    return net
end

for (k, v) in dsub100_to_dsub25
    net = get_net(:d100, k)
    if v == "GND"
        net.note = "GND"
        continue
    end
    net.d25_plug, net.d25_pin = parse.(Int, split(v, "-"))
end

for (k, v) in dsub100_to_electrode
    net = get_net(:d100, k)
    if v == "GND"
        @assert net.note == "GND"
        continue
    end
    net.electrode = v
end

for (k, v) in electrode_to_ao
    get_net(:electrode, k).dac = v
end

for i in 1:size(grid_to_electrode, 1)
    grid = grid_to_electrode[i, 1]
    electrode = grid_to_electrode[i, 2]
    if electrode == "SENSE1"
        get_net(:electrode, "TI1").grid = grid
        get_net(:electrode, "TV1").grid = grid
    elseif electrode == "SENSE2"
        get_net(:electrode, "TI2").grid = grid
        get_net(:electrode, "TV2").grid = grid
    elseif electrode == "GND"
        get_net(:grid, grid).note = "GND"
    else
        get_net(:electrode, electrode).grid = grid
    end
end

for i in 1:size(inter_pad_to_grid, 1)
    inter_pad = inter_pad_to_grid[i, 1]
    grid = inter_pad_to_grid[i, 2]
    net = get_net(:grid, grid)
    if net.electrode == "TI2" || net.electrode == "TV2"
        get_net(:electrode, "TI2").inter_pad = inter_pad
        get_net(:electrode, "TV2").inter_pad = inter_pad
    elseif net.electrode == "TI1" || net.electrode == "TV1"
        get_net(:electrode, "TI1").inter_pad = inter_pad
        get_net(:electrode, "TV1").inter_pad = inter_pad
    else
        net.inter_pad = inter_pad
    end
end

for d25_plug in 1:4
    for d25_pin in 1:25
        net = get_net(:d25_plug, d25_plug, :d25_pin, d25_pin)
        if isempty(net.d100)
            net.note = "GND"
        end
    end
end

sort!(all_nets, by=x->(x.electrode, x.grid, x.d100, x.note, x.dac, x.inter_pad))

function print_all(io::IO)
    println(io, ("D100 row,D100 pin in row,D100 pin,D25 plug,D25 pin,DAC chn," *
                 "Trap electrode,Package pad,CPGA grid,Note"))
    for net in all_nets
        print_csv(io, net)
        println(io)
    end
end

open(joinpath(@__DIR__, "../data/full_map.csv"), "w") do io
    print_all(io)
end
