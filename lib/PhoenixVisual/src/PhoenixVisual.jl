#!/usr/bin/julia

module PhoenixVisual

using EzXML
using Printf
using PyPlot
using PhoenixVoltages.Potentials

const svg_ns = ["svg" => "http://www.w3.org/2000/svg"]

mutable struct TrapTemplate
    const doc::EzXML.Document
    const svg::EzXML.Node
    const ionpos::EzXML.Node
    const y_ion::Float64
    const x_load::Float64
    const x_center::Float64
    title::EzXML.Node
    plotspace::EzXML.Node
    function TrapTemplate(doc::EzXML.Document)
        svg = root(doc)
        ionpos = findfirst("//svg:line[@id='IonPos']", svg, svg_ns)
        if ionpos === nothing
            error("Cannot find IonPos node")
        end
        y_ion = parse(Float64, ionpos["y1"])
        x_load = parse(Float64, ionpos["x2"])
        x_center = parse(Float64, ionpos["x1"])
        template = new(doc, svg, ionpos, y_ion, x_load, x_center)
        title = findfirst("//svg:text[@id='Title']", svg, svg_ns)
        if title !== nothing
            template.title = title
        end
        plotspace = findfirst("//svg:rect[@id='PlotSpace']", svg, svg_ns)
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
    return template.title
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

function num2str(v)
    v = round(v, digits=2)
    if isinteger(v)
        return string(Int(v))
    end
    return string(v)
end

function add_circle!(template, x, y, r)
    circ = addelement!(template.ionpos.parentnode, "circle")
    circ["cx"] = num2str(x_trap_to_svg(template, x))
    circ["cy"] = num2str(y_trap_to_svg(template, y))
    circ["r"] = num2str(scale_trap_to_svg(template, r))
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
    local last_y, pending
    is_pending = false
    for (x, y) in zip(xs, ys)
        x = x_trap_to_svg(template, x)
        if y > 1
            y = 1.0
        elseif y < 0
            y = 0.0
        end
        y = ytop + height * (1 - y)
        x = num2str(x)
        y = num2str(y)
        if @isdefined(last_y) && last_y == y
            is_pending = true
            pending = (x, y)
            continue
        end
        if is_pending
            push!(points, "$(pending[1]),$(pending[2])")
        end
        is_pending = false
        last_y = y
        push!(points, "$x,$y")
    end
    if is_pending
        push!(points, "$(pending[1]),$(pending[2])")
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

function val_to_color(val)
    @sprintf "%02x" clamp(round(Int, val * 0xff), UInt8)
end

function fill_electrodes!(template, values)
    electrodes = findall("//svg:path[@id]", template.svg, svg_ns)
    cmap = get_cmap("RdBu_r")
    for ele in electrodes
        id = split(ele["id"], "-")[1]
        value = get(values, id, 0)
        if value == 0
            continue
        end
        rgba = cmap(value / 2 + 0.5)
        color_str = "#" * val_to_color(rgba[1]) * val_to_color(rgba[2]) * val_to_color(rgba[3])
        ele["style"] = replace(ele["style"], r"fill:[^;]*"=>"fill:$color_str")
    end
end

function Base.write(io::IO, template::TrapTemplate)
    println(io, template.doc)
end

function _get_center(fits_cache, centers, xpos_um)
    xidx = Potentials.x_axis_to_index(fits_cache.solution, xpos_um ./ 1000)
    return (xidx, get(centers, xidx)...)
end

function render_frame(fits_cache::Potentials.FitCache, centers, electrodes, voltages;
                      finalize=true, # Set finalize to false if more objects needs to be added.
                      title=nothing,
                      # For position plotting
                      xpos_um=nothing,
                      pos_size=20, pos_fill="blueviolet",
                      # For potential plotting
                      plotx_ums=nothing, # Must be sorted
                      plot_yoffset=0.2, plot_yscale=0.2,
                      plot_axis=true, plot_axis_margin=200,
                      plot_axis_stroke="gainsboro",
                      plot_stroke="coral"
                      )
    template = get_template(title=title !== nothing, plot=plotx_ums !== nothing)

    if title !== nothing
        set_title!(template, title)
    end

    if xpos_um !== nothing
        rf_center_idx = _get_center(fits_cache, centers, xpos_um)
        ypos_um = Potentials.y_index_to_axis(fits_cache.solution,
                                             rf_center_idx[2]) * 1000
        c = add_circle!(template, xpos_um, ypos_um, pos_size)
        c["fill"] = pos_fill
    end

    electrode_voltages = zip(electrodes, voltages)

    if plotx_ums !== nothing
        function get_voltage(x_um)
            local pos = _get_center(fits_cache, centers, x_um)
            fit = Potentials.get_multi_electrodes(fits_cache, electrode_voltages,
                                                  (pos[3], pos[2], pos[1]))
            return fit[0, 0, 0]
        end
        ax_potential = get_voltage.(plotx_ums)

        if plot_axis
            line2 = PhoenixVisual.add_plotline!(template,
                                                [plotx_ums[1] - plot_axis_margin, plotx_ums[end] + plot_axis_margin],
                                                [plot_yoffset, plot_yoffset])
            line2["fill"] = "none"
            line2["stroke"] = plot_axis_stroke
        end

        line = PhoenixVisual.add_plotline!(template, plotx_ums,
                                           ax_potential .* plot_yscale .+ plot_yoffset)
        line["fill"] = "none"
        line["stroke"] = plot_stroke
    end

    voltage_map = Dict{String,Float64}()
    for (idx, v) in electrode_voltages
        v = v / 20
        for name in fits_cache.solution.electrode_names[idx]
            voltage_map[name] = v
        end
    end

    fill_electrodes!(template, voltage_map)

    if finalize
        finalize_svg!(template)
    end
    return template
end

end
