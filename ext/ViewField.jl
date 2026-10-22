# ViewField.jl
# A field on the plane of view_domain: scalars, edge and facet quantities, reconstructed
# with interpolate. Included from FITToolboxMakieExt.jl.

const _NodalEntity  = Union{FITToolbox.PrimalNode, FITToolbox.DualNode}
const _VectorEntity = Union{FITToolbox.PrimalEdge, FITToolbox.DualEdge,
                            FITToolbox.PrimalFacet, FITToolbox.DualFacet}

const _COMPLEX_MODES = ["instantaneous", "real part", "imaginary part", "magnitude", "phase"]

# --- Helpers without rendering (tested in test/viewer.jl) --------------------

# "name" => (data, entity) or "name" => (snapshots, entity, times) for each field,
# checked against the domain. A static field is a single snapshot without times.
function _field_specs(domain, fields)
    return map(collect(fields)) do f
        f isa Pair ||
            throw(ArgumentError("a field is given as \"name\" => (data, entity), got $(typeof(f))"))
        name, spec = f
        spec isa Tuple && length(spec) in (2, 3) ||
            throw(ArgumentError("field \"$name\": give (data, entity), e.g. (Φ, PrimalNode()), \
                                 or (snapshots, entity, times)"))
        data, entity = spec[1], spec[2]
        entity isa Union{_NodalEntity,_VectorEntity} ||
            throw(ArgumentError("field \"$name\": the entity must be PrimalNode(), DualNode(), \
                                 PrimalEdge(), DualEdge(), PrimalFacet() or DualFacet(), got $entity"))
        n = entity isa _NodalEntity ? domain.Np : 3 * domain.Np
        times = length(spec) == 3 ? collect(Float64, spec[3]) : nothing
        snapshots = times === nothing ? [data] : _snapshots(data)
        for s in snapshots
            length(s) == n && continue
            hint = times === nothing && eltype(data) <: AbstractVector ?
                   "; for snapshots give their times as well: (snapshots, entity, times)" : ""
            throw(ArgumentError("field \"$name\": $(nameof(typeof(entity))) data has length \
                                 $(length(s)), the domain needs $n$hint"))
        end
        times === nothing || length(times) == length(snapshots) ||
            throw(ArgumentError("field \"$name\": $(length(snapshots)) snapshots, but \
                                 $(length(times)) times"))
        (name = string(name), snapshots = snapshots, times = times, entity = entity,
         vector = entity isa _VectorEntity, complex = any(s -> eltype(s) <: Complex, snapshots))
    end
end

# The snapshots of a time series: a vector of vectors, or the columns of a matrix
_snapshots(data::AbstractMatrix) = [view(data, :, k) for k in axes(data, 2)]
_snapshots(data::AbstractVector{<:AbstractVector}) = collect(data)
_snapshots(data) = throw(ArgumentError("snapshots are a vector of vectors or a matrix with \
                                        one snapshot per column, got $(typeof(data))"))

# The snapshot closest to time t, and the time shown for a step
_nearest_step(times, t) = argmin(abs.(times .- t))
function _time_label(t)
    t == 0 && return "0 s"
    for (unit, name) in ((1.0, "s"), (1e-3, "ms"), (1e-6, "µs"), (1e-9, "ns"), (1e-12, "ps"))
        abs(t) >= unit && return "$(round(t / unit; sigdigits = 4)) $name"
    end
    return "$(round(t / 1e-15; sigdigits = 4)) fs"
end

# The quantities a field offers: its value, or the magnitude and components of a vector
_quantities(f) = f.vector ? ["magnitude", "x", "y", "z"] : ["value"]

# whether the shown values can be negative, i.e. need a range symmetric about zero
_signed(f, quantity, mode) = (f.complex && mode == "phase") ||       # ±180°, also for |vector|
    !(quantity == "magnitude" || (f.complex && mode == "magnitude"))

# The nodes of one axis with `refine - 1` points between each pair
function _refined(nodes, refine::Int)
    refine >= 1 || throw(ArgumentError("refine must be at least 1, got $refine"))
    refine == 1 && return collect(Float64, nodes)
    out = Float64[]
    for i in 1:length(nodes)-1, s in 0:refine-1
        push!(out, nodes[i] + s / refine * (nodes[i+1] - nodes[i]))
    end
    return push!(out, nodes[end])
end

# The two in-plane axes for a plane normal to `axis`, ascending, and their sample points
_plane_axes(axis::Int) = filter(!=(axis), (1, 2, 3))
function _plane_coords(domain, axis::Int, refine::Int)
    a, b = _plane_axes(axis)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    return _refined(nodes[a], refine), _refined(nodes[b], refine)
end

# The field at the sample points of the plane, from snapshot `step`: a matrix for
# scalars, an array with the three components along its third dimension for vectors
function _sample_plane(domain, f, axis::Int, pos, xs, ys, step::Int = 1)
    data = f.snapshots[clamp(step, 1, length(f.snapshots))]
    a = first(_plane_axes(axis))
    T = f.complex ? ComplexF64 : Float64
    out = f.vector ? Array{T}(undef, length(xs), length(ys), 3) :
                     Matrix{T}(undef, length(xs), length(ys))
    for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
        p = ntuple(d -> d == axis ? Float64(pos) : d == a ? x : y, 3)
        v = FITToolbox.interpolate(domain, f.entity, data, p...)
        f.vector ? (out[i, j, :] .= v) : (out[i, j] = v)
    end
    return out
end

# A complex value as shown: Re(z e^{jφ}) at phase φ (in degrees), or one of its parts
function _complex_part(z, mode, phase)
    mode == "instantaneous"  && return real(z * cis(deg2rad(phase)))
    mode == "real part"      && return real(z)
    mode == "imaginary part" && return imag(z)
    mode == "magnitude"      && return abs(z)
    mode == "phase"          && return rad2deg(angle(z))
    throw(ArgumentError("unknown mode \"$mode\""))
end

# The values to colour the plane with
function _display_values(samples, f, quantity, mode, phase)
    part(z) = f.complex ? _complex_part(z, mode, phase) : Float64(z)
    f.vector || return map(part, samples)
    if quantity == "magnitude"
        f.complex && mode == "phase" && return fill(NaN, size(samples)[1:2])   # a vector has no phase
        f.complex && mode == "magnitude" &&                  # the amplitude, over all phases
            return dropdims(sqrt.(sum(abs2, samples; dims = 3)); dims = 3)
        return dropdims(sqrt.(sum(abs2, part.(samples); dims = 3)); dims = 3)
    end
    return map(part, samples[:, :, findfirst(==(quantity), ("x", "y", "z"))])
end

# The largest value a point can show for this quantity and mode: the amplitude for the
# instantaneous value, so that the global range holds for every phase
function _bound(v, f, quantity, mode)
    if f.vector && quantity == "magnitude"
        f.complex || return sqrt(sum(abs2, v))
        mode in ("instantaneous", "magnitude") && return sqrt(sum(abs2, v))
        return sqrt(sum(abs2, mode == "real part" ? real.(v) : imag.(v)))
    end
    z = f.vector ? v[findfirst(==(quantity), ("x", "y", "z"))] : v
    f.complex || return abs(z)
    mode in ("instantaneous", "magnitude") && return abs(z)
    return abs(mode == "real part" ? real(z) : imag(z))
end

# The colour range from the largest value M (and the smallest positive one, for a log
# scale): symmetric about zero for signed values, from zero otherwise
function _range_from(M, m_pos, signed, log)
    (isfinite(M) && M > 0) || return signed ? (-1.0, 1.0) : (log ? (0.1, 1.0) : (0.0, 1.0))
    signed && return (-M, M)
    log && return (min(max(m_pos, M * 1e-6), M / 10), M)   # one to six decades below the maximum
    return (0.0, M)
end

# The colour range over the whole field and all its snapshots, from the field at the
# grid nodes (an evenly spaced subset of them above `max_points` in all)
function _global_range(domain, f, quantity, mode, log; max_points = 1_000_000)
    f.complex && mode == "phase" && return (-180.0, 180.0)
    s = max(1, ceil(Int, cbrt(domain.Np * length(f.snapshots) / max_points)))
    M, m_pos = 0.0, Inf
    for data in f.snapshots, k in 1:s:domain.Nw, j in 1:s:domain.Nv, i in 1:s:domain.Nu
        v = FITToolbox.interpolate(domain, f.entity, data,
                                   domain.nodes_u[i], domain.nodes_v[j], domain.nodes_w[k])
        b = _bound(v, f, quantity, mode)
        isfinite(b) || continue
        M = max(M, b)
        b > 0 && (m_pos = min(m_pos, b))
    end
    return _range_from(M, m_pos, _signed(f, quantity, mode), log)
end

# The colour range of the values on the plane only
function _slice_range(values, signed, log)
    vals = filter(isfinite, values)
    isempty(vals) && return _range_from(NaN, NaN, signed, log)
    pos = filter(>(0), vals)
    return _range_from(maximum(abs, vals), isempty(pos) ? NaN : minimum(pos), signed, log)
end

# The plane as a mesh through the sample points; vertex i + (j-1)·length(xs) is (xs[i], ys[j])
function _plane_mesh(axis::Int, pos, xs, ys)
    a = first(_plane_axes(axis))
    pts = [Point3f(ntuple(d -> d == axis ? pos : d == a ? x : y, 3)) for y in ys for x in xs]
    n = length(xs)
    idx(i, j) = i + (j - 1) * n
    faces = _TriFace[]
    for j in 1:length(ys)-1, i in 1:n-1
        push!(faces, _TriFace(idx(i, j), idx(i+1, j), idx(i+1, j+1)),
                     _TriFace(idx(i, j), idx(i+1, j+1), idx(i, j+1)))
    end
    return Makie.GeometryBasics.Mesh(pts, faces)
end

# --- The plane and its panel section -------------------------------------------

function _add_field!(ls, panel, row, domain, specs, axis, slider; refine = 1)
    any_complex = any(f -> f.complex, specs)
    any_series  = any(f -> f.times !== nothing, specs)

    Label(panel[row += 1, 1:3], "Field"; font = :bold, halign = :left)
    field_on = _checkbox(panel[row += 1, 1], false)
    field_menu = Menu(panel[row, 2:3]; options = [f.name for f in specs], default = specs[1].name)
    Label(panel[row += 1, 1], "show"; halign = :left)
    quantity_menu = Menu(panel[row, 2:3]; options = _quantities(specs[1]))
    mode_menu = phase = nothing
    if any_complex
        Label(panel[row += 1, 1], "complex"; halign = :left)
        mode_menu = Menu(panel[row, 2:3]; options = _COMPLEX_MODES, default = "instantaneous")
        phase = Slider(panel[row += 1, 2:3]; range = 0:5:360, startvalue = 0)
        Label(panel[row, 1], lift(p -> "$(p)°", phase.value); halign = :left)
    end
    time_slider = play = nothing
    if any_series                                    # the time itself is in the colour bar label
        first_series = specs[findfirst(f -> f.times !== nothing, specs)]
        play = _button(panel[row += 1, 1], "play")
        time_slider = Slider(panel[row, 2:3]; range = 1:max(2, length(first_series.times)),
                             startvalue = 1)
    end
    range_exp = Slider(panel[row += 1, 2:3]; range = -3:0.1:0, startvalue = 0)
    Label(panel[row, 1], lift(e -> "×$(round(10.0^e; sigdigits = 2))", range_exp.value);
          halign = :left)
    fit = _checkbox(panel[row += 1, 1], false)
    Label(panel[row, 2:3], "fit colours to plane"; halign = :left)
    logscale = _checkbox(panel[row += 1, 1], false)
    both = GridLayout(panel[row, 2:3]; halign = :left)
    Label(both[1, 1], "log scale"; halign = :left)
    viridis = _checkbox(both[1, 2], false)
    Label(both[1, 3], "viridis"; halign = :left)

    # untyped: the fields differ in type (real/complex data, entity), and an Observable
    # keeps the type of its first value
    field = Observable{Any}(specs[1])
    on(s -> field[] = specs[findfirst(f -> f.name == s, specs)], field_menu.selection)
    on(field) do f                                   # a vector or a scalar: other quantities
        quantity_menu.options[] = _quantities(f)
        quantity_menu.i_selected[] = 1
    end
    quantity = lift(q -> something(q, "value"), quantity_menu.selection)
    mode  = mode_menu === nothing ? Observable("instantaneous") : mode_menu.selection
    angle = phase === nothing ? Observable(0) : phase.value

    # time series: the slider steps through the snapshots of the field shown; switching
    # to another series keeps the time, at its nearest snapshot
    step = time_slider === nothing ? Observable(1) : time_slider.value
    playing = Observable(false)
    time_text = Observable("")                       # "t = 2.0 ns  (3/4)", or "static field"
    if time_slider !== nothing
        last_time = Ref(first_series.times[1])
        on(field) do f
            f.times === nothing && return
            time_slider.range[] = 1:max(2, length(f.times))
            set_close_to!(time_slider, _nearest_step(f.times, last_time[]))
        end
        onany(field, step) do f, k
            if f.times === nothing
                time_text[] = "static field"
            else
                k = clamp(k, 1, length(f.times))
                last_time[] = f.times[k]
                time_text[] = "t = $(_time_label(f.times[k]))  ($k/$(length(f.times)))"
            end
        end
        notify(step)
        on(p -> play.label[] = p ? "pause" : "play", playing)
        # steps on while playing, until the window is closed; a CairoMakie figure has no
        # window, so nothing runs in the tests
        on(play.clicks) do _
            playing[] = !playing[]
            playing[] || return
            @async try
                while playing[] && events(ls.scene).window_open[]
                    n = length(field[].snapshots)
                    n > 1 && set_close_to!(time_slider, mod1(time_slider.value[] + 1, n))
                    sleep(0.1)
                end
            catch err
                @error "playing the time series stopped" exception = (err, catch_backtrace())
            finally
                playing[] = false
            end
        end
    end

    # the samples, only while the field is shown; recomputed when the plane moves
    xs_ys = Observable(_plane_coords(domain, axis[], refine))
    samples = Observable{Any}(nothing)
    plane = Observable(_plane_mesh(axis[], slider.value[], xs_ys[]...))
    plane_at = Ref((axis[], slider.value[]))         # the mesh is rebuilt only when it moves
    onany(field_on.checked, field, axis, slider.value, step) do shown, f, a, pos, k
        shown || return
        xs, ys = _plane_coords(domain, a, refine)
        if (a, pos) != plane_at[]                   # not for a new snapshot
            xs_ys[] = (xs, ys)
            plane[] = _plane_mesh(a, pos, xs, ys)
            plane_at[] = (a, pos)
        end
        samples[] = (f, _sample_plane(domain, f, a, pos, xs, ys, k))   # with the field they belong to
    end

    values = Observable(fill(NaN32, prod(length, xs_ys[])))   # one per vertex, even hidden: GLMakie checks
    signed = Observable(true)
    use_log = Observable(false)
    ranges = Dict{Any,Any}()                         # global ranges, computed once each
    colorrange = Observable((-1.0, 1.0))
    onany(samples, quantity, mode, angle, fit.checked, logscale.checked,
          range_exp.value) do fs, q, m, φ, fitted, lg, e
        fs === nothing && return
        f, smp = fs
        q in _quantities(f) || return                # the quantity menu is being switched
        v = _display_values(smp, f, q, m, φ)
        sg = _signed(f, q, m)
        lg = lg && !sg                               # a log scale only for magnitudes
        rng = fitted ? _slice_range(v, sg, lg) :
              get!(() -> _global_range(domain, f, q, m, lg), ranges, (f.name, q, m, lg))
        # the range slider narrows it (the phase keeps its ±180°)
        f.complex && m == "phase" || (rng = (rng[1] * 10.0^e, rng[2] * 10.0^e))
        signed[] = sg
        # a log scale never sees a range from zero: the range goes first when switching
        # it on, last when switching it off
        lg ? (colorrange[] = rng; use_log[] = lg) : (use_log[] = lg; colorrange[] = rng)
        values[] = Float32.(vec(lg ? map(x -> x > 0 ? x : NaN, v) : v))
    end

    colormap = lift(v -> v ? :viridis : :jet, viridis.checked)
    scale = map!(l -> l ? log10 : identity, Observable{Function}(identity), use_log)   # two types
    plt = mesh!(ls, plane; color = values, colormap, colorrange, colorscale = scale,
                nan_color = RGBAf(0, 0, 0, 0), shading = Makie.NoShading,
                visible = field_on.checked)
    label = lift(field, quantity, mode, time_text) do f, q, m, t
        (f.vector ? "$(f.name): $q" * (f.complex ? ", $m" : "") :
                    f.name * (f.complex ? ": $m" : "")) * (f.times === nothing ? "" : ", $t")
    end
    Colorbar(panel[row += 1, 1:3]; vertical = false, flipaxis = false, colormap, colorrange,
             scale, label, height = 10)

    return row, (; on = field_on, field, quantity, mode, phase = angle, fit = fit.checked,
                 log = logscale.checked, viridis = viridis.checked, samples, values,
                 colorrange, colormap, signed, plot = plt, field_menu, quantity_menu, mode_menu,
                 phase_slider = phase, step, time_slider, time_text, play, playing,
                 range_slider = range_exp)
end
