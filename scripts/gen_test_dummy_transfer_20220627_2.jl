#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../lib"))

using PhoenixVoltages.OutputFiles

const mapfile = load_file(ARGS[1], MapFile)

const nelectrodes = length(mapfile.names)
const outputprefix = ARGS[2]

# To make the current loader happy, we need at least 4 lines
# LoadingJump, Loading, Center, Zeros.

const ntest_lines = 8
const lines = [fill(1.0, nelectrodes), fill(1.5, nelectrodes),
               fill(0.5, nelectrodes)]

function gen_test_line_n(n)
    n -= 1
    return [((i >> n) & 1) == 0 ? -0.1 : 0.1 for i in 0:(nelectrodes - 1)]
end
for i in 1:ntest_lines
    push!(lines, gen_test_line_n(i))
end
push!(lines, zeros(nelectrodes))

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
    for i in 1:ntest_lines
        start_name = i == 1 ? "Center" : "$(i - 1)um"
        println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"$(i + 1)\" startName=\"$(start_name)\" stopLine=\"$(i + 2)\" stopName=\"$(i)um\"/>")
    end
    println(io, "    <ShuttleEdge idleCount=\"500\" steps=\"5\" startLength=\"0\" stopLength=\"0\" startLine=\"$(ntest_lines + 2)\" startName=\"$(ntest_lines)um\" stopLine=\"$(ntest_lines + 3)\" stopName=\"Zeros\"/>")

    # XML footer
    print(io, """
  </ShuttlingGraph>
</VoltageAdjust>
""")
end
