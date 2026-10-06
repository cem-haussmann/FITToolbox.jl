# ViewDomain.jl
# Interactive 3D view of the objects in a domain. Included from FITToolboxMakieExt.jl.

# Above this many objects, the panel switches from one toggle per object to one
# toggle per object type; single objects are then hidden by clicking them.
const _MAX_OBJECT_TOGGLES = 15

# --- Helpers without rendering (tested in test/viewer.jl) --------------------

# 3D scenes need a depth buffer, which CairoMakie lacks: overlapping objects would
# be drawn in the wrong order.
function _check_3d_backend(backend = Makie.current_backend())
    backend === missing &&
        error("view_domain needs an active Makie backend with a depth buffer: \
               run `using GLMakie` (window) or `using WGLMakie` (notebook) first")
    nameof(backend) === :CairoMakie &&
        error("view_domain needs a backend with a depth buffer: CairoMakie cannot draw \
               overlapping 3D objects correctly. Run `using GLMakie` (window) or \
               `using WGLMakie` (notebook), or activate one with `GLMakie.activate!()`")
    return nothing
end

_type_name(obj) = string(nameof(typeof(obj)))

# Checkbox groups for the panel: (label, obj_ids, swatch colour or nothing).
# One group per object up to _MAX_OBJECT_TOGGLES objects, else one per object type.
function _toggle_groups(entries, colors)
    if length(entries) <= _MAX_OBJECT_TOGGLES
        return [(string(e.name, "  (", _type_name(e.object), ")"), [e.obj_id], colors[e.obj_id])
                for e in entries]
    end
    types = unique(_type_name(e.object) for e in entries)
    return map(types) do t
        ids = [e.obj_id for e in entries if _type_name(e.object) == t]
        (string(t, " (", length(ids), ")"), ids, nothing)
    end
end

# The colour given to the object, or the next one from Makie's palette
function _object_colors(entries)
    palette = Makie.wong_colors()
    colors = Dict{Int,Any}()
    for (k, e) in enumerate(entries)
        m = e.object isa FITToolbox.AbstractSolid ? e.object.material : nothing
        colors[e.obj_id] = m !== nothing && m.color !== nothing ?
                           FITToolbox._color_hex(m.color) : palette[mod1(k, length(palette))]
    end
    return colors
end

_axis_index(::DirX) = 1
_axis_index(::DirY) = 2
_axis_index(::DirZ) = 3

# The positions the cut can take along axis 1, 2 or 3: the grid nodes
_cut_positions(domain, axis::Int) = (domain.nodes_u, domain.nodes_v, domain.nodes_w)[axis]

# Clip plane that hides everything beyond `pos` along `axis`. Makie clips a point
# when dot(normal, p) < distance, so the normal points towards the part that is kept.
function _clip_plane(axis::Int, pos)
    n = Vec3f(ntuple(i -> i == axis ? -1f0 : 0f0, 3))
    return Plane3f(Point3f(ntuple(i -> i == axis ? Float32(pos) : 0f0, 3)), n)
end

# Outline of the cut plane: a closed rectangle across the domain at `pos`
function _cut_rectangle(lo, hi, axis::Int, pos)
    a, b = filter(!=(axis), 1:3)
    corner(s, t) = Point3f(ntuple(i -> i == axis ? pos : i == a ? (s ? hi[a] : lo[a]) :
                                                   (t ? hi[b] : lo[b]), 3))
    return [corner(false, false), corner(true, false), corner(true, true),
            corner(false, true), corner(false, false)]
end

# The loop as a polygon in its plane
function _loop_points(loop::FITToolbox.CircularLoop; n = 128)
    a = _axis_index(loop.normal)
    b, c = filter(!=(a), 1:3)
    return map(range(0, 2π; length = n)) do t
        p = collect(Float32.(loop.center))
        p[b] += loop.radius * cos(t)
        p[c] += loop.radius * sin(t)
        Point3f(p)
    end
end

# Grid lines in the cut plane at `pos` along `axis`: one line per node of each of the
# two in-plane axes, as pairs of points for linesegments
function _cut_grid_segments(domain, axis::Int, pos)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    a, b = filter(!=(axis), 1:3)
    pt(x, y) = Point3f(ntuple(i -> i == axis ? pos : i == a ? x : y, 3))
    segs = Point3f[]
    for x in nodes[a]
        push!(segs, pt(x, nodes[b][1]), pt(x, nodes[b][end]))
    end
    for y in nodes[b]
        push!(segs, pt(nodes[a][1], y), pt(nodes[a][end], y))
    end
    return segs
end

# The primal edges a loop source occupies, from the same discretisation get_source
# uses: the edges as pairs of points, their centres, and the direction of the current
# on each (edge vector times the ±1 orientation, independent of the amplitude)
function _loop_edges(domain, loop::FITToolbox.CircularLoop)
    J = FITToolbox._circular_loop_unit_source(domain, loop, loop.normal)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    edges = (domain.edges_u, domain.edges_v, domain.edges_w)
    segs, mids, dirs = Point3f[], Point3f[], Vec3f[]
    for p in findall(!iszero, J)
        c, q = divrem(p - 1, domain.Np)          # component and node, both 0-based
        k, r = divrem(q, domain.Nu * domain.Nv)
        j, i = divrem(r, domain.Nu)
        ijk  = (i + 1, j + 1, k + 1)
        from = Point3f(ntuple(d -> nodes[d][ijk[d]], 3))
        step = Vec3f(ntuple(d -> d == c + 1 ? edges[d][ijk[d]] : 0.0, 3))
        push!(segs, from, from + step)
        push!(mids, from + step / 2)
        push!(dirs, sign(J[p]) * step)
    end
    return segs, mids, dirs
end

# --- Voxel view ------------------------------------------------------------

# The object each cell belongs to: the obj_id of the solid that wrote it last, 0 for
# the background. The history is replayed into a scratch domain on the same grid with
# every solid's σ set to its obj_id, so overlaps resolve exactly as in the material.
# Values written directly into domain.material belong to no object and are not shown.
function _owner_labels(domain)
    scratch = FITToolbox.create_domain(domain.edges_u, domain.edges_v, domain.edges_w; units = "m")
    # the objects warned about clipping etc. when they were created
    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        for e in FITToolbox.list_objects(domain).entries
            e.object isa FITToolbox.AbstractSolid || continue
            FITToolbox._apply!(scratch, _with_material(e.object, FITToolbox.Material(e.obj_id, 1.0, 1.0)))
        end
    end
    return round.(Int, view(scratch.material, :, :, :, 1))
end

# a copy of a solid with another material
_with_material(o, m) = typeof(o)((f === :material ? m : getfield(o, f) for f in fieldnames(typeof(o)))...)

# The six faces of a cell: corners wound counter-clockwise seen from outside, so the
# normals point out (corner c = 1 + dx + 2dy + 4dz), and the step to the neighbour
const _CELL_FACES = (((1, 5, 7, 3), (-1, 0, 0)), ((2, 4, 8, 6), (1, 0, 0)),
                     ((1, 2, 6, 5), (0, -1, 0)), ((3, 7, 8, 4), (0, 1, 0)),
                     ((1, 3, 4, 2), (0, 0, -1)), ((5, 6, 8, 7), (0, 0, 1)))

const _TriFace = Makie.GeometryBasics.GLTriangleFace

# Stands in for an object without cells (e.g. all of it beyond the cut): one
# degenerate triangle, the same mesh type as the real ones
const _EMPTY_VOXELS = Makie.GeometryBasics.Mesh(fill(Point3f(0), 3), [_TriFace(1, 2, 3)])

# The cells to keep for a cut at `pos` along `axis`: those that end at or before it.
# The slider stops on nodes, so the cut always falls on a cell boundary.
function _cut_cells(domain, axis::Int, pos)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)[axis]
    return (axis, searchsortedlast(nodes, pos + 1e-9 * (nodes[end] - nodes[1])) - 1)
end

"""
Surface meshes of the labelled cells, one per label. A face is drawn only where the
neighbour carries another label, so the meshes grow with the surface of the objects,
not with their volume. The corners come from the node coordinates, so graded grids
are drawn as they are. With `keep = (axis, m)`, cells beyond index `m` along `axis`
count as empty: a cut then shows a filled cross-section instead of a hollow shell.
"""
_voxel_meshes(domain, label; keep = nothing) =
    keep === nothing ? _voxel_meshes(domain, label, 0, 0) : _voxel_meshes(domain, label, keep...)

function _voxel_meshes(domain, label::AbstractArray{Int,3}, keep_axis::Int, keep_m::Int)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    n = size(label)
    lab(i, j, k) = (1 <= i <= n[1] && 1 <= j <= n[2] && 1 <= k <= n[3] &&
                    (keep_axis == 0 || (i, j, k)[keep_axis] <= keep_m)) ? label[i, j, k] : 0
    pts = Dict{Int,Vector{Point3f}}()
    fcs = Dict{Int,Vector{_TriFace}}()
    @inbounds for k in 1:n[3], j in 1:n[2], i in 1:n[1]
        id = lab(i, j, k)
        id == 0 && continue
        corner(c) = Point3f(nodes[1][i + ((c - 1) & 1)], nodes[2][j + ((c - 1) >> 1 & 1)],
                            nodes[3][k + ((c - 1) >> 2 & 1)])
        for (quad, (di, dj, dk)) in _CELL_FACES
            lab(i + di, j + dj, k + dk) == id && continue
            P = get!(Vector{Point3f}, pts, id)
            F = get!(Vector{_TriFace}, fcs, id)
            b = length(P)
            for c in quad
                push!(P, corner(c))
            end
            push!(F, _TriFace(b + 1, b + 2, b + 3), _TriFace(b + 1, b + 3, b + 4))
        end
    end
    return Dict(id => Makie.GeometryBasics.Mesh(pts[id], fcs[id]) for id in keys(pts))
end

# --- Panel widgets ---------------------------------------------------------

# A checkbox framed and, when ticked, filled in `color`: the object's own colour in the
# object list, the accent colour for settings
_checkbox(pos, checked; color = Makie.COLOR_ACCENT[]) =
    Checkbox(pos; checked, size = 18, checkboxstrokewidth = 2, checkboxcolor_checked = color,
             checkboxstrokecolor_checked = color, checkboxstrokecolor_unchecked = color)

_button(pos, label) = Button(pos; label, cornerradius = 10, buttoncolor = RGBf(0.90, 0.92, 0.96))

# Back to the camera at the start: Makie's default viewing direction, centred on the scene
function _reset_view!(scene)
    cam = Makie.cameracontrols(scene)
    Makie.update_cam!(scene, cam, Vec3f(3), Vec3f(0), Vec3f(0, 0, 1))
    Makie.center!(scene)
    return nothing
end

# --- Drawing ---------------------------------------------------------------

_draw_object!(scene, b::FITToolbox.Brick, color; kw...) =
    mesh!(scene, Rect3f(Point3f(b.origin), Vec3f(b.lengths)); color, kw...)

_draw_object!(scene, c::FITToolbox.Cylinder, color; kw...) =
    mesh!(scene, Makie.GeometryBasics.normal_mesh(Makie.Tessellation(
                     Makie.GeometryBasics.Cylinder(Point3f(c.base), Point3f(_cylinder_top(c)),
                                                   Float32(c.radius)), 64));
          color, kw...)

_cylinder_top(c) = ntuple(d -> d == _axis_index(c.axis) ? c.base[d] + c.height : c.base[d], 3)

_draw_object!(scene, s::FITToolbox.Sphere, color; kw...) =
    mesh!(scene, Makie.GeometryBasics.normal_mesh(
                     Makie.Tessellation(Makie.Sphere(Point3f(s.center), Float32(s.radius)), 64));
          color, kw...)

# the caller sets the linewidth: it changes when the loop is shown on its grid edges
_draw_object!(scene, c::FITToolbox.CircularLoop, color; kw...) =
    lines!(scene, _loop_points(c); color, kw...)

# Walk up from a picked primitive to the plot drawn for an object
function _picked_id(plot_ids, plt)
    while plt !== nothing
        haskey(plot_ids, plt) && return plot_ids[plt]
        plt = plt.parent isa Makie.AbstractPlot ? plt.parent : nothing
    end
    return nothing
end

function FITToolbox.view_domain(domain::FITToolbox.FITDomain; size = (1100, 700))
    _check_3d_backend()
    fig, _ = _view_figure(domain; size)
    return fig
end

# The figure with all its controls, and the state they act on, for the tests:
# CairoMakie cannot draw it, but can build it and run every callback.
function _view_figure(domain::FITToolbox.FITDomain; size = (1100, 700))

    entries = FITToolbox.list_objects(domain).entries
    colors  = _object_colors(entries)
    lo = (domain.nodes_u[1], domain.nodes_v[1], domain.nodes_w[1])
    hi = (domain.nodes_u[end], domain.nodes_v[end], domain.nodes_w[end])

    fig = Figure(; size)
    ls  = LScene(fig[1, 1]; show_axis = true)

    # the domain, never clipped, so the cut keeps its frame
    wireframe!(ls, Rect3f(Point3f(lo), Vec3f(hi .- lo)); color = :gray60)

    # one visibility switch per object; toggles, clicks and "Show all" all set these
    shown  = Dict(e.obj_id => Observable(true) for e in entries)
    planes = Observable(Plane3f[])
    plot_ids = IdDict{Any,Int}()
    # switched by the "View" menu and the "Loops" toggles in the panel
    voxel_mode = Observable(false)
    edges_on   = Observable(false)
    arrows_on  = Observable(false)
    # loops show the edges they occupy when asked to, and always in the voxel view,
    # which is the discretised view as a whole
    edges_shown = lift(|, edges_on, voxel_mode)
    for e in entries
        id, c = e.obj_id, colors[e.obj_id]
        if e.object isa FITToolbox.CircularLoop
            # the exact circle, thin while the edges are shown on top of it
            plt = _draw_object!(ls, e.object, c; visible = shown[id], clip_planes = planes,
                                linewidth = lift(on -> on ? 1 : 3, edges_shown))
            plot_ids[plt] = id
            segs, mids, dirs = _loop_edges(domain, e.object)
            isempty(segs) && continue
            plt = linesegments!(ls, segs; color = c, linewidth = 4, clip_planes = planes,
                                visible = lift(&, shown[id], edges_shown))
            plot_ids[plt] = id
            # arrows along the edges, 80 % of an edge long, sized to the smallest
            # cell so they stay visible next to the edge lines on any grid
            h = minimum(d -> maximum(abs, d), dirs)
            plt = arrows3d!(ls, mids, 0.8f0 .* dirs; color = c, align = :center,
                            markerscale = h, tiplength = 0.45, tipradius = 0.25,
                            shaftradius = 0.08, minshaftlength = 0.1, clip_planes = planes,
                            visible = lift(&, shown[id], edges_shown, arrows_on))
            plot_ids[plt] = id
        else
            plt = _draw_object!(ls, e.object, c; clip_planes = planes,
                                visible = lift((s, v) -> s && !v, shown[id], voxel_mode))
            plot_ids[plt] = id
        end
    end

    # the voxels of each solid; built when the voxel view is first opened and rebuilt
    # when the cut changes (cells beyond the cut count as empty, so the cut face is
    # drawn), hence no clip planes here
    voxels = Dict(e.obj_id => Observable(_EMPTY_VOXELS)
                  for e in entries if e.object isa FITToolbox.AbstractSolid)
    for (id, m) in voxels
        plt = mesh!(ls, m; color = colors[id], visible = lift(&, shown[id], voxel_mode))
        plot_ids[plt] = id
    end

    # --- panel ---
    panel = GridLayout(fig[1, 2]; tellheight = false, valign = :top)
    row = 0
    Label(panel[row += 1, 1], "View"; font = :bold, halign = :left)
    view_menu = Menu(panel[row, 2:3]; options = ["objects", "voxels"], default = "objects")
    connect!(voxel_mode, lift(==("voxels"), view_menu.selection))
    Label(panel[row += 1, 1:3], "Objects"; font = :bold, halign = :left)
    isempty(entries) && Label(panel[row += 1, 1:3], "none"; halign = :left)

    toggles = Checkbox[]
    for (label, ids, swatch) in _toggle_groups(entries, colors)
        t = _checkbox(panel[row += 1, 1], true; color = something(swatch, :gray50))
        Label(panel[row, 2:3], label; halign = :left)
        on(checked -> foreach(id -> shown[id][] = checked, ids), t.checked)
        push!(toggles, t)
    end

    info = Observable(isempty(entries) ? " " : "Click an object to hide it")
    Label(panel[row += 1, 1:3], info; halign = :left, tellwidth = false)
    buttons  = GridLayout(panel[row += 1, 1:3])
    show_all = _button(buttons[1, 1], "Show all")
    reset    = _button(buttons[1, 2], "Reset view")
    on(_ -> _reset_view!(ls.scene), reset.clicks)
    on(show_all.clicks) do _
        foreach(t -> t.checked[] = true, toggles)
        foreach(o -> o[] = true, values(shown))
        info[] = "Click an object to hide it"
    end

    # --- cut ---
    Label(panel[row += 1, 1:3], "Cut"; font = :bold, halign = :left)
    cut_on = _checkbox(panel[row += 1, 1], false)
    axis_menu = Menu(panel[row, 2:3]; options = ["x", "y", "z"], default = "z")
    pos0 = _cut_positions(domain, 3)
    slider = Slider(panel[row += 1, 1:3]; range = pos0, startvalue = pos0[cld(length(pos0), 2)])
    Label(panel[row += 1, 1:3], lift(v -> "position: $(round(v; sigdigits = 4)) m", slider.value);
          halign = :left)

    axis = lift(s -> findfirst(==(s), ["x", "y", "z"]), axis_menu.selection)
    on(axis) do a
        pos = _cut_positions(domain, a)
        slider.range[] = pos
        set_close_to!(slider, pos[cld(length(pos), 2)])
    end
    onany(cut_on.checked, axis, slider.value) do active, a, pos
        planes[] = active ? [_clip_plane(a, pos)] : Plane3f[]
    end

    labels = Ref{Union{Nothing,Array{Int,3}}}(nothing)     # computed on first use
    function update_voxels!()
        voxel_mode[] || return
        labels[] === nothing && (labels[] = _owner_labels(domain))
        keep = cut_on.checked[] ? _cut_cells(domain, axis[], slider.value[]) : nothing
        meshes = _voxel_meshes(domain, labels[]; keep)
        for (id, m) in voxels
            m[] = get(meshes, id, _EMPTY_VOXELS)
        end
    end
    onany((_...) -> update_voxels!(), voxel_mode, cut_on.checked, axis, slider.value)

    # the grid in the cut plane, shown while the cut is active
    grid_on = _checkbox(panel[row += 1, 1], false)
    Label(panel[row, 2:3], "grid in cut plane"; halign = :left)
    linesegments!(ls, lift((a, pos) -> _cut_grid_segments(domain, a, pos), axis, slider.value);
                  color = (:gray20, 0.6), linewidth = 0.75,
                  visible = lift(&, cut_on.checked, grid_on.checked))

    # --- loops on the grid edges they occupy ---
    if any(e -> e.object isa FITToolbox.CircularLoop, entries)
        Label(panel[row += 1, 1:3], "Loops"; font = :bold, halign = :left)
        t = _checkbox(panel[row += 1, 1], false)
        Label(panel[row, 2:3], "on grid edges"; halign = :left)
        connect!(edges_on, t.checked)
        t = _checkbox(panel[row += 1, 1], false)
        Label(panel[row, 2:3], "current direction"; halign = :left)
        connect!(arrows_on, t.checked)
    end
    lines!(ls, lift((a, pos) -> _cut_rectangle(lo, hi, a, pos), axis, slider.value);
           color = :black, linestyle = :dash, visible = cut_on.checked)

    # --- click to hide: a click is a press and release without moving, so
    # dragging still rotates the camera ---
    press_pos = Ref(Vec2f(NaN))
    on(events(fig).mousebutton) do event
        event.button == Mouse.left || return Consume(false)
        mouse = events(fig).mouseposition[]
        if event.action == Mouse.press
            press_pos[] = mouse
        elseif event.action == Mouse.release && is_mouseinside(ls.scene) &&
               hypot((mouse .- press_pos[])...) < 3
            plt, _ = pick(ls.scene)
            id = _picked_id(plot_ids, plt)
            if id !== nothing
                shown[id][] = false
                e = only(x for x in entries if x.obj_id == id)
                info[] = "hidden: $(e.name) (obj_id $id)"
            end
        end
        return Consume(false)
    end

    colsize!(fig.layout, 2, Fixed(280))
    return fig, (; ls, shown, voxels, voxel_mode, planes, plot_ids, info)
end
