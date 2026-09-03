# Plotting.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 06/08/2026
# Collapsed from three near-identical methods into one, 11/08/2026

using CairoMakie

# For a cut with normal `n`, return everything that differs between the three
# orientations: the two in-plane node vectors, the two in-plane edge-centre
# vectors, which components of the field lie in the plane, and the axis labels.
_slice_axes(::X, c) = (c.nodes_v, c.nodes_w, c.edges_v_center, c.edges_w_center, 2, 3, "y", "z")
_slice_axes(::Y, c) = (c.nodes_u, c.nodes_w, c.edges_u_center, c.edges_w_center, 1, 3, "x", "z")
_slice_axes(::Z, c) = (c.nodes_u, c.nodes_v, c.edges_u_center, c.edges_v_center, 1, 2, "x", "y")

# Rebuild a full (x, y, z) triple from the two in-plane coordinates and `pos`.
_at(::X, pos, a, b) = (pos, a, b)
_at(::Y, pos, a, b) = (a, pos, b)
_at(::Z, pos, a, b) = (a, b, pos)

# Index of the entry of the sorted vector `v` nearest to `x`. Binary search
# rather than argmin(abs.(v .- x)), which allocates on every streamline step.
@inline function _nearest(v, x)
    i = searchsortedfirst(v, x)
    i <= 1 && return 1
    i > length(v) && return length(v)
    return (x - v[i-1]) < (v[i] - x) ? i - 1 : i
end


function plot_nodal_values(config, ::Primal, data::AbstractVector, normal::Normals;
                           pos = 0, units = "m", logscale = false,
                           label = "Potential Φ [V]", clip_range = nothing,
                           plot_negative_gradient = false,
                           density = 1.2, stepsize = 0.5, maxsteps = 1000,
                           gridsize = (32, 32))

    length(data) == config.Np ||
        throw(DimensionMismatch("data has length $(length(data)), expected Np = $(config.Np)"))

    pos = Float64(pos) * check_units(units)
    an, bn, ac, bc, c1, c2, alab, blab = _slice_axes(normal, config)

    # Sample the potential at every node of the cut plane.
    slice = Matrix{Float32}(undef, length(an), length(bn))
    for (jb, b) in enumerate(bn), (ia, a) in enumerate(an)
        slice[ia, jb] = interpolate(config, PrimalNode(), data, _at(normal, pos, a, b)...; units = "m")
    end

    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1];
              xlabel = "$alab-direction (m)",
              ylabel = "$blab-direction (m)",
              aspect = DataAspect())

    # On a log colour scale the auto-range must ignore zero and negative values,
    # which a potential routinely has.
    c_range = if !isnothing(clip_range)
        clip_range
    elseif logscale
        pos_vals = filter(>(0), slice)
        isempty(pos_vals) && throw(ArgumentError("logscale = true but no positive values in the slice"))
        (minimum(pos_vals), maximum(pos_vals))
    else
        (minimum(slice), maximum(slice))
    end

    hm = heatmap!(ax, collect(an), collect(bn), slice;
                  colormap = cgrad(:jet, 256),
                  interpolate = true,
                  colorrange = c_range,
                  (logscale ? (; colorscale = log10) : (;))...)

    if plot_negative_gradient
        field = -(G(config) * data)

        # Field components live on primal edges, so sample at edge centres.
        vec_field = Array{Float32,3}(undef, length(ac), length(bc), 3)
        for (jb, b) in enumerate(bc), (ia, a) in enumerate(ac)
            e = interpolate(config, PrimalEdge(), field, _at(normal, pos, a, b)...; units = "m")
            vec_field[ia, jb, 1] = e[1]
            vec_field[ia, jb, 2] = e[2]
            vec_field[ia, jb, 3] = e[3]
        end

        as, bs = collect(ac), collect(bc)

        # c1 and c2 pick the two components tangent to the cut plane; the third
        # is discarded, so off-centre slices show a projection.
        f(a, b) = Vec2f(vec_field[_nearest(as, a), _nearest(bs, b), c1],
                        vec_field[_nearest(as, a), _nearest(bs, b), c2])

        streamplot!(ax, f, extrema(as), extrema(bs);
                    density, stepsize, maxsteps, gridsize,
                    color = _ -> RGBAf(1, 1, 1, 1),
                    linewidth = 1.2)
    end

    Colorbar(fig[1, 2], hm; label = label)
    return fig
end