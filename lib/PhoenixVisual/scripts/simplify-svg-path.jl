#!/usr/bin/julia

using EzXML
using Parsers

const svg_ns = ["svg" => "http://www.w3.org/2000/svg"]
const svg = readxml(ARGS[1])
const svgroot = root(svg)

function read_command(io)
    while !eof(io)
        c = read(io, Char)
        isspace(c) || return c
    end
    return '\0'
end

function read_numbers!(io, num_buf)
    empty!(num_buf)
    while !eof(io)
        c = peek(io, Char)
        if isspace(c) || c == ','
            read(io, Char)
            continue
        end
        if isletter(c)
            read(io, Char)
            return c
        end
        push!(num_buf, Parsers.parse(Float64, io))
    end
    return '\0'
end

function parse_path(d)
    points = Tuple{Float64,Float64}[]
    io = IOBuffer(d)
    next_cmd = read_command(io)
    num_buf = Float64[]
    cur_pos = (0.0, 0.0)
    while next_cmd != '\0'
        cmd = next_cmd
        next_cmd = read_numbers!(io, num_buf)
        nargs = length(num_buf)
        if cmd == 'z' || cmd == 'Z'
            @assert nargs == 0
            @assert next_cmd == '\0'
            break
        elseif cmd == 'm'
            @assert nargs % 2 == 0
            for i in 1:2:nargs
                cur_pos = cur_pos .+ (num_buf[i], num_buf[i + 1])
                push!(points, cur_pos)
            end
        elseif cmd == 'M'
            @assert nargs % 2 == 0
            for i in 1:2:nargs
                cur_pos = (num_buf[i], num_buf[i + 1])
                push!(points, cur_pos)
            end
        elseif cmd == 'h'
            for i in 1:nargs
                cur_pos = (cur_pos[1] + num_buf[i], cur_pos[2])
                push!(points, cur_pos)
            end
        elseif cmd == 'H'
            for i in 1:nargs
                cur_pos = (num_buf[i], cur_pos[2])
                push!(points, cur_pos)
            end
        elseif cmd == 'v'
            for i in 1:nargs
                cur_pos = (cur_pos[1], cur_pos[2] + num_buf[i])
                push!(points, cur_pos)
            end
        elseif cmd == 'V'
            for i in 1:nargs
                cur_pos = (cur_pos[1], num_buf[i])
                push!(points, cur_pos)
            end
        elseif cmd == 'l'
            @assert nargs % 2 == 0
            for i in 1:2:nargs
                cur_pos = cur_pos .+ (num_buf[i], num_buf[i + 1])
                push!(points, cur_pos)
            end
        elseif cmd == 'L'
            @assert nargs % 2 == 0
            for i in 1:2:nargs
                cur_pos = (num_buf[i], num_buf[i + 1])
                push!(points, cur_pos)
            end
        else
            error("Unknown command $cmd in $d")
        end
    end
    return points
end

"""
Make the longest segment the last one
"""
function cycle_path!(points)
    maxdist = 0.0
    maxidx = 0
    for i in 1:length(points)
        p0 = i == 1 ? points[end] : points[i - 1]
        p1 = points[i]
        dist = hypot((p0 .- p1)...)
        if dist > maxdist
            maxdist = dist
            maxidx = i
        end
    end
    @assert maxidx != 0
    # maxidx will become the fist point
    points .= [@view(points[maxidx:end]); @view(points[1:maxidx - 1])]
    return
end

function collapse_closeby_points!(points, threshold)
    p0 = prev_p = points[1]
    prev_i = 1
    vec0 = prev_vec = prev_p .- points[end]
    for i in 2:length(points)
        p = points[i]
        dist = hypot((p .- prev_p)...)
        if dist <= threshold
            continue
        end
        # Now we need to decide where to put the previous points
        vec = p .- points[i - 1]
        # The previous line was prev_p + prev_vec * t1
        # The next line is p + vec * t2
        # dot(ovec, vec) == 0
        ovec = (vec[2], -vec[1])
        t1 = sum((p .- prev_p) .* ovec) / sum(prev_vec .* ovec)
        p_cross = prev_p .+ prev_vec .* t1
        dist = hypot((p_cross .- prev_p)...)
        if dist <= threshold * 2
            points[prev_i] = p_cross
        end
        prev_p = p
        prev_i += 1
        prev_vec = vec
        points[prev_i] = p
    end
    # check the position of the last point
    ovec0 = (vec0[2], -vec0[1])
    t1 = sum((p0 .- prev_p) .* ovec0) / sum(prev_vec .* ovec0)
    p_cross = prev_p .+ prev_vec .* t1
    dist = hypot((p_cross .- prev_p)...)
    if dist <= threshold * 2
        points[prev_i] = p_cross
    end
    if length(points) != prev_i
        resize!(points, prev_i)
    end
    return
end

function round_points!(points)
    for i in 1:length(points)
        p = points[i]
        points[i] = (round(p[1], digits=3), round(p[2], digits=3))
    end
end

function simplify_path!(points, threshold)
    cycle_path!(points)
    collapse_closeby_points!(points, threshold)
    round_points!(points)
end

function num2str(v)
    v = round(v, digits=3)
    if isinteger(v)
        return string(Int(v))
    end
    return string(v)
end

function encode_path(points)
    io = IOBuffer()
    prev_p = points[1]
    print(io, "m$(num2str(prev_p[1])),$(num2str(prev_p[2]))")
    prev_cmd = 'l'
    for i in 2:length(points)
        p = points[i]
        if p[1] == prev_p[1]
            abs_s = num2str(p[2])
            rel_s = num2str(p[2] - prev_p[2])
            if length(abs_s) >= length(rel_s)
                print(io, "v$(rel_s)")
            else
                print(io, "V$(abs_s)")
            end
        elseif p[2] == prev_p[2]
            abs_s = num2str(p[1])
            rel_s = num2str(p[1] - prev_p[1])
            if length(abs_s) >= length(rel_s)
                print(io, "h$(rel_s)")
            else
                print(io, "H$(abs_s)")
            end
        else
            abs_s = num2str(p[1]) * "," * num2str(p[2])
            rel_s = num2str(p[1] - prev_p[1]) * "," * num2str(p[2] - prev_p[2])
            if length(abs_s) >= length(rel_s)
                print(io, "l$(rel_s)")
            else
                print(io, "L$(abs_s)")
            end
        end
        prev_p = p
    end
    print(io, 'z')
    return String(take!(io))
end

id = 0

for p in findall("//svg:path", svgroot, svg_ns)
    points = parse_path(p["d"])
    if startswith(p["id"], "path")
        global id
        simplify_path!(points, 0.05)
        p["id"] = "path$(id)"
        id += 1
    else
        simplify_path!(points, 0.1)
    end
    p["d"] = encode_path(points)
end
open(ARGS[1], "w") do io
    print(io, svgroot)
end
