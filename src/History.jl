#History.jl was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 05/10/2026

Base.propertynames(d::FITDomain, private::Bool=false) =
    private ? fieldnames(typeof(d)) : filter(!=(:_history), fieldnames(typeof(d)))

function _record_create!(domain::FITDomain, object::AbstractFITObject; name::String)
    h = domain._history
    h.obj_counter += 1
    push!(h.steps, CreateObject(h.obj_counter, name, object))
    return h.obj_counter
end

function _record_delete!(domain::FITDomain, obj_id::Int)
    h = domain._history
    obj_id in _active_ids(h) ||
        throw(ArgumentError("no object with obj_id $obj_id in the domain"))
    push!(h.steps, DeleteObject(obj_id))
    return nothing
end

function _active_ids(h::DomainHistory)
    ids = Int[]
    for s in h.steps
        s isa CreateObject && push!(ids, s.obj_id)
        s isa DeleteObject && filter!(!=(s.obj_id), ids)
    end
    return ids
end

function _active_entries(h::DomainHistory)
    active = Set(_active_ids(h))
    return CreateObject[s for s in h.steps if s isa CreateObject && s.obj_id in active]
end

const _DEFAULT_NAME = r"^obj\d+$"

# final name for the next object: "obj<id>" if none is given; user names must be unique
function _resolve_name(domain::FITDomain, name::AbstractString)
    h = domain._history
    isempty(name) && return "obj$(h.obj_counter + 1)"
    occursin(_DEFAULT_NAME, name) &&
        throw(ArgumentError("names of the form \"obj<number>\" are reserved for unnamed objects, got \"$name\""))
    any(e -> e.name == name, _active_entries(h)) &&
        throw(ArgumentError("an object named \"$name\" already exists in the domain"))
    return String(name)
end

# the existing object with this name or id
function _entry(domain::FITDomain, name::AbstractString)
    for e in _active_entries(domain._history)
        e.name == name && return e
    end
    throw(ArgumentError("no object named \"$name\" in the domain"))
end

function _entry(domain::FITDomain, obj_id::Int)
    for e in _active_entries(domain._history)
        e.obj_id == obj_id && return e
    end
    throw(ArgumentError("no object with obj_id $obj_id in the domain"))
end

# Rebuild domain.material from the history: background first, then every existing
# solid in creation order, so later objects overwrite earlier ones as before.
# Values written directly into domain.material are lost; warn=true says so. Callers
# pass warn = !_material_matches_history(domain), checked before the history changes.
function _replay!(domain::FITDomain; warn::Bool=false)
    warn && @warn "Rebuilding the material from the history: values written directly into \
           domain.material are discarded. Use create_brick!, create_sphere!, … to keep them."
    md = domain.material
    md[:, :, :, 1] .= domain.σ
    md[:, :, :, 2] .= domain.ε_r
    md[:, :, :, 3] .= domain.μ_r

    active = Set(_active_ids(domain._history))
    for step in domain._history.steps
        step isa CreateObject            || continue   # skip DeleteObject steps
        step.obj_id in active            || continue   # skip objects deleted later
        step.object isa AbstractSolid    || continue   # sources do not paint material
        _apply!(domain, step.object)
    end
    return domain
end

# Throws if a source does not fit the grid of `domain`, before it is recorded or a grid
# with it is handed out. By building its source vector; edge sources and charges have
# cheaper checks of their own, without the vector of the whole grid.
_check_source(domain::FITDomain, src::AbstractSource) = (_discretize_source(domain, src); nothing)

# Every source of the history must be discretisable on the grid of `domain`. A loop that
# no longer fits within nodes[1] … nodes[end-2], or a polyline whose points all snap to
# one node, throws here, where the grid changes or is read from a file, and not later in
# get_source or view_domain.
function _check_sources(domain::FITDomain)
    for e in _active_entries(domain._history)
        e.object isa AbstractSource || continue
        try
            _check_source(domain, e.object)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("source \"$(e.name)\" (obj_id $(e.obj_id)) does not fit this \
                                 grid: $(err.msg)"))
        end
    end
    return nothing
end

"""
    remove_object!(domain, obj_id)

Remove the object `obj_id` from `domain` and rebuild the material from the history.
The removal is recorded as a step of its own, so `undo!` brings the object back.
Throws an `ArgumentError` if no object with that id exists.
"""
function remove_object!(domain::FITDomain, obj_id::Int)
    _entry(domain, obj_id)              # throws for unknown or already deleted ids
    clean = _material_matches_history(domain)
    _record_delete!(domain, obj_id)
    _replay!(domain; warn=!clean)
    return nothing
end

remove_object!(domain::FITDomain, name::AbstractString) =
    remove_object!(domain, _entry(domain, name).obj_id)

"""
    undo!(domain)

Revert the last step in the history of `domain`: a created object is removed again,
a removed object comes back. The material is rebuilt from the history. Object ids
are never handed out twice, so the next object still gets a new `obj_id`.

Throws an `ArgumentError` if the history is empty.
"""
function undo!(domain::FITDomain)
    h = domain._history
    isempty(h.steps) &&
        throw(ArgumentError("nothing to undo: the history of this domain is empty"))
    clean = _material_matches_history(domain)
    pop!(h.steps)
    _replay!(domain; warn=!clean)
    return nothing
end
# true if domain.material is exactly what the history produces,
# i.e. nothing was written into it directly
function _material_matches_history(domain::FITDomain)
    tmp = deepcopy(domain)
    # an internal check: objects must not warn again, e.g. about a clipped sphere
    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        _replay!(tmp)
    end
    return tmp.material == domain.material
end

"""
    change_resolution(domain, resolution; units="m") -> FITDomain
    change_resolution(domain, Edges_U, Edges_V, Edges_W; units="m") -> FITDomain

A new domain with the same extent and background material as `domain`, and the
full history of `domain` replayed on it.

The first form gives an equidistant grid: `resolution` is one value for all
directions or three values (u, v, w), each of which must divide the extent of
`domain` in its direction. The second form takes the cell widths per
direction, as `create_domain(Edges_U, Edges_V, Edges_W)` does, for non-equidistant
grids; each must add up to the extent of `domain` in that direction. All lengths
are given in `units`.

`domain` itself is not changed, and both domains have independent histories from
then on. Values written directly into `domain.material` are not part of the
history; if there are any, a warning says they are not carried over.

Throws an `ArgumentError` if a resolution does not divide the extent of `domain` or
the edges do not add up to it, or if an object cannot be created on the new grid, e.g.
a brick thinner than one cell, or a loop source that no longer fits within
`nodes[1] … nodes[end-2]`.
Source vectors obtained earlier with `get_source` are not updated; call
`get_source` on the new domain.
"""
function change_resolution(domain::FITDomain, Edges_U, Edges_V, Edges_W; units="m")
    unitToMeter = check_units(units)
    for (dir, edges, nodes) in (("u", Edges_U, domain.nodes_u),
                                ("v", Edges_V, domain.nodes_v),
                                ("w", Edges_W, domain.nodes_w))
        L = nodes[end] - nodes[1]
        isapprox(sum(edges) * unitToMeter, L; rtol=1e-9) ||
            throw(ArgumentError("the $dir edges add up to $(sum(edges) * unitToMeter) m, \
                                 but the domain is $L m long in $dir"))
    end
    _material_matches_history(domain) ||
        @warn "domain.material contains values that were not created through create_* \
               functions; they are not carried over to the new grid"

    new = create_domain(Edges_U, Edges_V, Edges_W; units=units,
                        σ=domain.σ, ε_r=domain.ε_r, μ_r=domain.μ_r)
    # copy the history: steps are immutable, so sharing them is safe
    append!(new._history.steps, domain._history.steps)
    new._history.obj_counter = domain._history.obj_counter
    _replay!(new; warn=false)
    _check_sources(new)
    return new
end

function change_resolution(domain::FITDomain, resolution; units="m")
    resolution isa Real && (resolution = (resolution, resolution, resolution))
    length(resolution) == 3 ||
        throw(ArgumentError("resolution needs 1 or 3 values, got $(length(resolution))"))
    unitToMeter = check_units(units)
    extents = (domain.nodes_u[end] - domain.nodes_u[1],
               domain.nodes_v[end] - domain.nodes_v[1],
               domain.nodes_w[end] - domain.nodes_w[1])
    # The extent stays, so each resolution must divide it: one error that says so, not a
    # warning per direction from _cell_count and then an error about the edges
    function cells(dir, L, r)
        h = r * unitToMeter
        h > 0 || throw(ArgumentError("resolution must be positive, got $r $units in $dir"))
        n = round(Int, L / h)
        n >= 1 && isapprox(n * h, L; rtol=1e-9) ||
            throw(ArgumentError("resolution $r $units does not divide the extent of the domain \
                                 in $dir, $(_short(L / unitToMeter)) $units \
                                 ($(round(L / h; sigdigits=3)) cells): give one that does, or \
                                 the cell widths, change_resolution(domain, Edges_U, Edges_V, Edges_W)"))
        return fill(h, n)
    end
    edges = [cells(dir, L, r) for (dir, L, r) in zip(("u", "v", "w"), extents, resolution)]
    return change_resolution(domain, edges...; units="m")
end

"""
    get_source(domain, name)   -> Vector
    get_source(domain, obj_id) -> Vector

Source vector of one source in `domain`, discretized on the current grid: for a loop
from `create_circular_loop_source!` or `create_rectangular_loop_source!` and a wire from
`create_polyline_source!` the current vector, length `3Np` in the edge layout, scaled by
its `current`; for a point charge from `create_electric_charge!` the charge
vector, length `Np`, with the charge in the dual volume of the nearest node. The element
type is `ComplexF64` for a complex current, else `Float64`.

Call it again after `change_resolution`, `remove_object!` or `undo!`: a vector
obtained earlier is not updated by these.

Throws an `ArgumentError` if no such object exists or it is not a source.
"""
get_source(domain::FITDomain, name::AbstractString) = get_source(domain, _entry(domain, name).obj_id)

# Every current on edges (circular and rectangular loops, polylines) is kept sparse: each
# type provides _unit_edges, the edges it occupies with their current per unit, a
# Dict edge => ±1 (±2 where a polyline runs along an edge twice), checked to be closed.
# Creating, checking, measuring and drawing a source use these alone, so their work and
# memory grow with the length of the wire, not with the grid. The source vector is
# written from them once, with the current: the only vector of the whole grid.
function _discretize_source(config::FITDomain, src::AbstractEdgeSource)
    J = zeros(typeof(src.current), 3 * config.Np)
    for (edge, value) in _unit_edges(config, src)
        J[edge] = src.current * value
    end
    return J
end

_check_source(config::FITDomain, src::AbstractEdgeSource) = (_unit_edges(config, src); nothing)

function get_source(domain::FITDomain, obj_id::Int)
    e = _entry(domain, obj_id)
    e.object isa AbstractSource ||
        throw(ArgumentError("object $obj_id \"$(e.name)\" is a $(nameof(typeof(e.object))), not a source"))
    return _discretize_source(domain, e.object)
end

# Result of list_objects: behaves like a Vector{CreateObject}, prints as a table
struct ObjectList <: AbstractVector{CreateObject}
    entries::Vector{CreateObject}
end

Base.size(l::ObjectList) = size(l.entries)
Base.getindex(l::ObjectList, i::Int) = l.entries[i]

"""
    list_objects(domain) -> ObjectList

All objects that currently exist in `domain`, in creation order; removed objects
are left out. Displays as a table. Index or iterate it like a vector; each entry
has the fields `obj_id`, `name` and `object`.

# Example

```julia
list_objects(domain)                     # show the table
remove_object!(domain, "coil")           # by name
J = get_source(domain, "coil")           # one source by name
```
"""
list_objects(domain::FITDomain) = ObjectList(_active_entries(domain._history))

function Base.show(io::IO, ::MIME"text/plain", l::ObjectList)
    n = length(l)
    print(io, n, n == 1 ? " object" : " objects")
    n == 0 && return
    println(io, ":")
    header = ("obj_id", "type", "name", "geometry (m)", "material")
    rows = [(string(e.obj_id), string(nameof(typeof(e.object))),
             isempty(e.name) ? "-" : e.name,
             _describe_geometry(e.object), _describe_material(e.object)) for e in l]
    w = [maximum(textwidth, (header[c], (r[c] for r in rows)...)) for c in 1:5]
    line(r) = join((rpad(r[c], w[c]) for c in 1:5), "  ")
    print(io, "  ", rstrip(line(header)))
    for r in rows
        print(io, "\n  ", rstrip(line(r)))
    end
end

# for display only: 12 significant digits hide float noise like 0.5000000000000003
_short(x::Real)    = round(x; sigdigits=12)
_short(x::Complex) = complex(_short(real(x)), _short(imag(x)))
_short(t::Tuple)   = map(_short, t)

_describe_material(m::Material) =
    "σ=$(_short(m.σ)) ε_r=$(_short(m.ε_r)) μ_r=$(_short(m.μ_r))" * (m.color === nothing ? "" : " color=$(m.color)")
