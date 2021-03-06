#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages.Outputs

const orig_comp_file = ARGS[1]
const trans_file = ARGS[2]
const new_comp_file = ARGS[3]

const comp = load_file(orig_comp_file, CompensationFile)
const trans = load_file(trans_file, TransferFile)
comp.term_values[9] = trans.line_values[177]
comp.term_values[10] = trans.line_values[207]
write_file(new_comp_file, comp)
