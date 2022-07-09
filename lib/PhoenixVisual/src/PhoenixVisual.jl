#!/usr/bin/julia

module PhoenixVisual

using EzXML

mutable struct TrapTemplate
    const doc::EzXML.Document
    const ionpos::EzXML.Node
    const y_ion::Float64
    const x_load::Float64
    const x_center::Float64
    title::EzXML.Node
    plotspace::EzXML.Node
    function TrapTemplate(doc::EzXML.Document)
        svg = root(doc)
        ns = ["svg" => "http://www.w3.org/2000/svg"]
        ionpos = findfirst("//svg:line[@id='IonPos']", svg, ns)
        if ionpos === nothing
            error("Cannot find IonPos node")
        end
        y_ion = parse(Float64, ionpos["y1"])
        x_load = parse(Float64, ionpos["x2"])
        x_center = parse(Float64, ionpos["x1"])
        template = new(doc, ionpos, y_ion, x_load, x_center)
        title = findfirst("//svg:text[@id='Title']", svg, ns)
        if title !== nothing
            template.title = title
        end
        plotspace = findfirst("//svg:rect[@id='PlotSpace']", svg, ns)
        if plotspace !== nothing
            template.plotspace = plotspace
        end
        return template
    end
end

function get_template(; title=false, plot=false)
    if plot
        name = "template-extra-space.svg"
    elseif title
        name = "template-title.svg"
    else
        name = "template-raw.svg"
    end
    return TrapTemplate(readxml(joinpath(@__DIR__, "../data", name)))
end

function set_title!(template, title)
    if !isdefined(template, :title)
        error("SVG template does not support title")
    end
    template.title.content = title
    return
end

function scale_trap_to_svg(template, v)
    return v * (template.x_center - template.x_load) / 3045
end

function x_trap_to_svg(template, x)
    return scale_trap_to_svg(template, x) + template.x_center
end

function y_trap_to_svg(template, y)
    return scale_trap_to_svg(template, y) + template.y_ion
end

function add_circle!(template, x, y, r)
    circ = addelement!(template.ionpos.parentnode, "circle")
    circ["cx"] = x_trap_to_svg(template, x)
    circ["cy"] = y_trap_to_svg(template, y)
    circ["r"] = scale_trap_to_svg(template, r)
    return circ
end

function add_plotline!(template, xs, ys)
    if !isdefined(template, :plotspace)
        error("SVG template does not support plotting")
    end
    plotspace = template.plotspace
    line = addelement!(plotspace.parentnode, "polyline")
    ytop = parse(Float64, plotspace["y"])
    height = parse(Float64, plotspace["height"])
    points = String[]
    for (x, y) in zip(xs, ys)
        x = x_trap_to_svg(template, x)
        if y > 1
            y = 1.0
        elseif y < 0
            y = 0.0
        end
        y = ytop + height * (1 - y)
        push!(points, "$x,$y")
    end
    line["points"] = join(points, " ")
    return line
end

function finalize_svg!(template)
    unlink!(template.ionpos)
    if isdefined(template, :plotspace)
        unlink!(template.plotspace)
    end
end

function Base.write(io::IO, template::TrapTemplate)
    println(io, template.doc)
end

end
