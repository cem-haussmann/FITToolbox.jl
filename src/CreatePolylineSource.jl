#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 21/10/2026

# The points are a tuple, not a vector: like every object in the history the polyline
# must compare by value (two equal histories are ==), which tuples do and vectors do not.
struct PolylineSource{T<:Union{Float64,ComplexF64}, N} <: AbstractEdgeSource
    points::NTuple{N,NTuple{3,Float64}}   # m, as given; the nodes are chosen when discretised
    closed::Bool                           # true: the last point joins back to the first
    current::T                             # A; complex for phasors in the frequency domain
end

# Int/Float32 → Float64, any Complex → ComplexF64
PolylineSource(points, closed, current::Real) =
    PolylineSource{Float64, length(points)}(Tuple(points), closed, current)
PolylineSource(points, closed, current::Complex) =
    PolylineSource{ComplexF64, length(points)}(Tuple(points), closed, current)

_describe_geometry(w::PolylineSource) =
    "$(length(w.points)) points from $(_short(first(w.points))) to $(_short(last(w.points)))" *
    (w.closed ? ", closed" : "")
_describe_material(w::PolylineSource) = "I=$(_short(w.current)) A"

"""
    create_polyline_source!(domain, points; closed = false, units = "m", current = 1.0,
                            name = "") -> obj_id

Add a wire carrying `current` (in A) along a polyline to `domain` and return its
`obj_id`; its discrete current vector comes from `get_source(domain, obj_id)` or
`get_source(domain, name)`, length `3Np` in the edge layout. `points` is a vector of
points `(u, v, w)` in `units`; the wire runs straight from each to the next, and with
`closed = true` from the last back to the first. A straight wire is two points.

On the grid, each point snaps to its nearest node, and each straight piece becomes a
path of primal edges between them: walking along the line, the path always steps along
the axis whose next cell-centre plane the line crosses first (on a tie the lower axis,
u before v before w). The path is therefore connected by construction, one edge per
step, and runs through the nodes nearest to the line, within half a cell diagonal of
it, on graded grids too. Along a grid axis it is exact; a diagonal piece becomes a
staircase. Every edge carries ±`current` in the direction of travel.

The wire is divergence free everywhere except at its ends: an open polyline has the
divergence `+current` at its first node, where the current enters the grid, and
`-current` at its last. For static problems such a source has no solution unless it is
closed (`closed = true`, or ending where it started) or ends on a boundary that takes the
current, such as PEC; in the time domain an open wire is a valid antenna. Creating an
open polyline therefore prints an info message naming the two end nodes.

`current` may be complex; the source vector is then `ComplexF64`.

# Errors

Throws an `ArgumentError` for fewer than two points (three if `closed`), for a point
outside the domain, if all points snap to the same node, or if `name` is already used.
Nothing is recorded then.

# Example

```julia
# a straight wire across the domain, from (0.1, 0.2, 0.5) m to (0.9, 0.7, 0.5) m
create_polyline_source!(domain, [(0.1, 0.2, 0.5), (0.9, 0.7, 0.5)]; current = 10.0,
                        name = "wire")
# a closed triangle in the plane z = 0.3 m
create_polyline_source!(domain, [(0.2, 0.2, 0.3), (0.8, 0.2, 0.3), (0.5, 0.7, 0.3)];
                        closed = true, name = "triangle")
```
"""
function create_polyline_source!(domain::FITDomain, points::AbstractVector; closed::Bool=false,
                                 units="m", current=1.0, name::String="")
    unitToMeter = check_units(units)
    all(p -> length(p) == 3, points) ||
        throw(ArgumentError("every point needs 3 coordinates (u, v, w)"))
    name = _resolve_name(domain, name)
    wire = PolylineSource([ntuple(d -> Float64(p[d]) * unitToMeter, 3) for p in points],
                          closed, current)
    _discretize_source(domain, wire)     # throws on invalid input, before anything is recorded
    first_node, last_node = _polyline_nodes(domain, wire)[[1, end]]
    if !closed && first_node != last_node
        at(n) = _short(ntuple(d -> (domain.nodes_u, domain.nodes_v, domain.nodes_w)[d][n[d]], 3))
        @info "Open polyline source: the current enters the grid at node $(at(first_node)) m \
               and leaves it at $(at(last_node)) m, so it is not divergence free there. For \
               static problems close it (closed = true) or end it on a boundary that takes \
               the current."
    end
    return _record_create!(domain, wire; name)
end

"""
    get_source_ends(domain, name)   -> (start = p, stop = q) or nothing
    get_source_ends(domain, obj_id) -> (start = p, stop = q) or nothing

The two end nodes of a polyline source from `create_polyline_source!`, as node indices
in `1:Np`, the layout of `get_index_entity(domain, PrimalNode(), …)` and of nodal
vectors such as a potential. `start` is the node where the current enters the grid,
`stop` where it leaves: the nodes the first and the last point snap to, independent of
the sign of the current. For a positive current the divergence of its source vector is
`+current` at `start` and `-current` at `stop`.

Returns `nothing` if the polyline has no ends: if it is closed, or ends where it started.
Throws an `ArgumentError` if no such object exists or it is not a polyline source.

# Example

```julia
ends = get_source_ends(domain, "wire")
get_position_of_index(domain, PrimalNode(), ends.start)     # its coordinates
```
"""
get_source_ends(domain::FITDomain, name::AbstractString) =
    get_source_ends(domain, _entry(domain, name).obj_id)

function get_source_ends(domain::FITDomain, obj_id::Int)
    e = _entry(domain, obj_id)
    e.object isa PolylineSource ||
        throw(ArgumentError("object $obj_id \"$(e.name)\" is a $(nameof(typeof(e.object))), \
                             not a polyline source"))
    path = _polyline_nodes(domain, e.object)
    path[1] == path[end] && return nothing
    node(n) = n[1] + (n[2]-1)*domain.Nu + (n[3]-1)*domain.Nu*domain.Nv
    return (start = node(path[1]), stop = node(path[end]))
end

# The nodes the points snap to, in order; for a closed polyline the first again at the end
function _polyline_nodes(config::FITDomain, w::PolylineSource)
    nodes = (config.nodes_u, config.nodes_v, config.nodes_w)
    n = length(w.points)
    n >= (w.closed ? 3 : 2) ||
        throw(ArgumentError("a$(w.closed ? " closed" : "") polyline needs at least \
                             $(w.closed ? 3 : 2) points, got $n"))
    snapped = map(w.points) do p
        all(d -> _within(p[d], nodes[d]), 1:3) ||
            throw(ArgumentError("polyline point $p m lies outside the domain"))
        ntuple(d -> _find_index(nodes[d], p[d]), 3)
    end
    return w.closed ? [collect(snapped); [snapped[1]]] : collect(snapped)
end

function _unit_source(config::FITDomain, w::PolylineSource)
    path = _polyline_nodes(config, w)
    all(==(path[1]), path) &&
        throw(ArgumentError("the polyline collapses on the grid: all its points snap to \
                             the same node"))
    J = zeros(Float64, 3 * config.Np)
    for k in 1:length(path)-1
        _walk_edges!(J, config, path[k], path[k+1])
    end

    # Divergence free except at the ends of an open polyline: S̃J = +1 where the current
    # enters, −1 where it leaves. True by construction; checked as for the loops.
    expected = zeros(config.Np)
    node_index(n) = n[1] + (n[2]-1)*config.Nu + (n[3]-1)*config.Nu*config.Nv
    if path[1] != path[end]
        expected[node_index(path[1])]   =  1.0
        expected[node_index(path[end])] = -1.0
    end
    residual = norm(get_divergence(config, Dual()) * J - expected)
    residual > 1e-12 * norm(J) &&
        @warn "The polyline source is not connected: its divergence is wrong away from the ends" residual
    return J
end

# Add the path of edges from node P to node Q to J, ±1 in the direction of travel. The
# path follows the straight line between the two nodes: at every step it moves one edge
# along the axis whose next cell-centre plane the line crosses first, so it always runs
# through the nodes nearest to the line. On a tie the lower axis goes first.
function _walk_edges!(J, config::FITDomain, P, Q)
    nodes   = (config.nodes_u, config.nodes_v, config.nodes_w)
    centres = (config.edges_u_center, config.edges_v_center, config.edges_w_center)
    Nu, Nv, Np = config.Nu, config.Nv, config.Np
    xa = ntuple(d -> nodes[d][P[d]], 3)
    xb = ntuple(d -> nodes[d][Q[d]], 3)
    cur = collect(P)
    while Tuple(cur) != Tuple(Q)
        axis, t_first = 0, Inf
        for d in 1:3
            cur[d] == Q[d] && continue
            # the cell-centre plane between the current node and the next one along d
            c = Q[d] > cur[d] ? centres[d][cur[d]] : centres[d][cur[d]-1]
            t = (c - xa[d]) / (xb[d] - xa[d])
            t < t_first && ((axis, t_first) = (d, t))     # strict: ties keep the lower axis
        end
        step = Q[axis] > cur[axis] ? 1 : -1
        lower = step > 0 ? cur[axis] : cur[axis] - 1      # the edge starts at its lower node
        ijk = ntuple(d -> d == axis ? lower : cur[d], 3)
        J[(axis-1)*Np + ijk[1] + (ijk[2]-1)*Nu + (ijk[3]-1)*Nu*Nv] += step
        cur[axis] += step
    end
    return J
end
