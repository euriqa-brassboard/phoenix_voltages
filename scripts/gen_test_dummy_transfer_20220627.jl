#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages.OutputFiles

const mapfile = load_file(ARGS[1], MapFile)

const nelectrodes = length(mapfile.names)
const outputprefix = ARGS[2]

# To make the current loader happy, we need at least 4 lines
# LoadingJump, Loading, Center, Zeros.
# We'll set their voltages to 1.0, 1.5, 0.5, 0 respectively.
# (This voltage will be applied to all of the electrodes.)

const lines = [fill(1.0, nelectrodes), fill(1.5, nelectrodes),
               fill(0.5, nelectrodes), zeros(nelectrodes)]

const transfer_file = TransferFile(mapfile, lines)
write_file("$(outputprefix).txt", transfer_file)
open("$(outputprefix).xml", "w") do io
    # XML header
    print(io, """
<?xml version="1.0" ?>
<VoltageAdjust>
  <ShuttlingGraph>
""")

    println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"0\" startName=\"LoadingJump\" stopLine=\"1\" stopName=\"Loading\"/>")
    println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"1\" startName=\"Loading\" stopLine=\"2\" stopName=\"Center\"/>")
    println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"2\" startName=\"Center\" stopLine=\"3\" stopName=\"Zeros\"/>")

    # XML footer
    print(io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
end
