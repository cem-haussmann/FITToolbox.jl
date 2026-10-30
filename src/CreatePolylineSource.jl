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
`get_source(domain, name)`, length `3Np` in the edge layout. `points` is a vector (or a
tuple) of points `(u, v, w)` in `units`; the wire runs straight from each to the next, and with
`closed = true` from the last back to the first. A straight wire is two points.

On the grid, each point snaps to its nearest node, and each straight piece becomes a
path of primal edges between them: walking along the line, the path always steps along
the axis whose next cell-centre plane the line crosses first (on a tie the lower axis,
u before v before w). The path is therefore connected by construction, one edge per
step, and runs through the nodes nearest to the line, within half a cell diagonal of
it, on graded grids too. Along a grid axis it is exact; a diagonal piece becomes a
staircase. Every edge carries ±`current` in the direction of travel; `get_source_path`
gives the edges in their order along the wire.

The wire is divergence free everywhere except at its ends: an open polyline has the
divergence `+current` at its first node, where the current enters the grid, and
`-current` at its last. For static problems such a source has no solution unless it is
closed (`closed = true`, or ending where it started) or ends on a boundary that takes the
current, such as PEC; in the time domain an open wire is a valid antenna. Creating an
open polyline therefore prints an info message naming the two end nodes.

`current` may be complex; the source vector is then `ComplexF64`.

# Errors

Throws an `ArgumentError` for fewer than two points (three if `closed`), for a point
outside the domain, if all points snap to the same node, if the polyline cancels itself
(every edge walked as often back as forth, e.g. there and back along one line), if
`current` is not finite, or if `name` is already used. Nothing is recorded then.

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
function create_polyline_source!(domain::FITDomain, points::Union{AbstractVector,Tuple};
                                 closed::Bool=false,
                                 units="m", current=1.0, name::AbstractString="")
    unitToMeter = check_units(units)
    all(p -> length(p) == 3, points) ||
        throw(ArgumentError("every point needs 3 coordinates (u, v, w)"))
    name = _resolve_name(domain, name)
    wire = PolylineSource([ntuple(d -> Float64(p[d]) * unitToMeter, 3) for p in points],
                          closed, current)
    _check_source(domain, wire)          # throws on invalid input, before anything is recorded
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
    get_source_ends(domain, name)   -> (start, stop, start_node, stop_node) or nothing
    get_source_ends(domain, obj_id) -> (start, stop, start_node, stop_node) or nothing

The two ends of a polyline source from `create_polyline_source!`: `start` is where the
current enters the grid, `stop` where it leaves, as positions `(u, v, w)` in m, like the
ends of a helix from `get_helix_ends`. They are the nodes the first and the last point
snap to, independent of the sign of the current, so `start_node` and `stop_node` give
them as node indices in `1:Np` as well, the layout of
`get_index_entity(domain, PrimalNode(), …)` and of nodal vectors such as a potential.
For a positive current the divergence of the source vector is `+current` at
`start_node` and `-current` at `stop_node`.

Returns `nothing` if the polyline has no ends: if it is closed, or ends where it started.
Throws an `ArgumentError` if no such object exists or it is not a polyline source. The
whole path between the ends, edge by edge in order, comes from `get_source_path`.

# Example

```julia
ends = get_source_ends(domain, "wire")
ends.start                      # (u, v, w) in m
Φ[ends.start_node]              # the potential there
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
    at(n) = (domain.nodes_u[n[1]], domain.nodes_v[n[2]], domain.nodes_w[n[3]])
    return (start = at(path[1]), stop = at(path[end]),
            start_node = node(path[1]), stop_node = node(path[end]))
end

"""
    get_source_path(domain, name)   -> (edges, signs, lengths, nodes, pieces, planned_lengths)
    get_source_path(domain, obj_id) -> (edges, signs, lengths, nodes, pieces, planned_lengths)

The edges of a polyline source from `create_polyline_source!` in the order the current
runs along them, e.g. to build a transmission-line model of the wire on the grid:

- `edges`: their indices in the edge layout, `1:3Np` (see `get_curl`);
- `signs`: `+1.0` where the current runs along the edge's axis, `-1.0` against it;
- `lengths`: their lengths in m, which differ between the axes and on a graded grid;
- `nodes`: the nodes the path passes, in `1:Np`, one more than there are edges: edge
  `edges[m]` runs from `nodes[m]` to `nodes[m+1]`. An open polyline runs from the
  `start_node` to the `stop_node` of `get_source_ends`; a closed one ends where it starts;
- `pieces`: the straight piece of the polyline each edge belongs to, `k` for the piece
  from point `k` to point `k+1`; the closing piece of a closed polyline is the last;
- `planned_lengths`: the length of the planned wire each edge stands for, in m (below).

The path is the walk of `get_source` without its steps back along the edge just walked.
At a sharp corner the next piece of a polyline can run back over the end of the previous
one: those steps cancel in the source vector, and in a path they would be a dead end,
there and back. So every edge of the source vector is in the path once, and the edges
with their signs add up to it exactly, `Σₘ signs[m] e[edges[m]] = get_source(…) / current`.
A closed polyline whose last step goes back along its first edge starts at the next node.

The path is a staircase, longer than the line it follows (see `get_source_length`): a
straight piece of direction `Δ` becomes `‖Δ‖₁/‖Δ‖₂` times as long, `√2` at 45° in a
plane, `√3` along a space diagonal. `planned_lengths` undo this piece by piece: the length
of each piece, from the points as given, is shared among the edges walked for it in
proportion to their lengths, so that they add up to it exactly. A piece whose points snap
to the same node walks no edge; its length goes to the next piece that does, so a curve
drawn with points closer than a cell keeps its length too. Its single edges then get
uneven shares, which even out over a few cells; drawn with points about a cell apart,
each piece's edges share its length evenly. A per-length quantity of the
real wire, e.g. the inductance and capacitance of a transmission-line model, belongs on
the edges with `planned_lengths`: the wire then keeps its electrical length piece by
piece, also where it is partly straight and partly diagonal.

Where parts of the wire lie within a cell of each other, e.g. a leg running back next to
another, or a corner sharper than the grid can follow, their edges cancel, and their
length goes with them: no other part of the wire is stretched to make up for it.
`sum(planned_lengths)` then falls short of the planned length of `get_source_length` by
just that, and a warning says so when it is more than 5 %. A finer grid there, or the
parts at least a cell apart, keep them.

Unlike `get_source`, it builds no vector of the whole grid: its work and memory grow with
the length of the path, not with the grid, as do those of `create_polyline_source!` and
`get_source_length` for a polyline.

Throws an `ArgumentError` if no such object exists or it is not a polyline source, if it
does not fit the grid (as `get_source` does), or if it runs along an edge twice on the
grid, overlapping itself, since its edges then have no single order.

# Example

```julia
p = get_source_path(domain, "dipole")
M = length(p.edges)
Γ = sparse(1:M, p.edges, p.signs, M, 3 * domain.Np)    # segment m to its edge, along the path
Γ' * ones(M) == get_source(domain, "dipole")            # for a current of 1 A
L = L′ .* p.planned_lengths                             # per segment, as on the real wire
points = [get_position_of_index(domain, PrimalNode(), q)[1:3] for q in p.nodes]
```
"""
get_source_path(domain::FITDomain, name::AbstractString) =
    get_source_path(domain, _entry(domain, name).obj_id)

function get_source_path(domain::FITDomain, obj_id::Int)
    e = _entry(domain, obj_id)
    e.object isa PolylineSource ||
        throw(ArgumentError("object $obj_id \"$(e.name)\" is a $(nameof(typeof(e.object))), \
                             not a polyline source"))
    walk = _polyline_walk(domain, e.object)            # the checks of get_source
    allunique(walk.edges) ||
        throw(ArgumentError("the polyline \"$(e.name)\" runs along an edge twice on the grid, \
                             overlapping itself, so its edges have no single order"))
    # Only steps that cancel are left out, so the path holds the edges of the source vector
    # with their signs. True by construction; anything else is a bug in the package.
    length(walk.net) == length(walk.edges) &&
        all(get(walk.net, edge, 0.0) == sign for (edge, sign) in zip(walk.edges, walk.signs)) ||
        error("the path of the polyline \"$(e.name)\" differs from its source vector: this is \
               a bug in FITToolbox, please report it")
    lengths = [_edge_length(domain, edge) for edge in walk.edges]
    planned_lengths, planned = _planned_lengths(domain, e.object, walk.pieces, lengths)
    # what the grid cannot hold: parts of the wire whose edges cancel
    lost = 1 - sum(planned_lengths) / planned
    lost > 0.05 &&
        @warn "polyline \"$(e.name)\": $(round(Int, 100 * lost)) % of its planned length \
               cancels on the grid, where parts of it lie within a cell of each other; \
               planned_lengths cover only the rest"
    return (edges = walk.edges, signs = walk.signs, lengths = lengths, nodes = walk.nodes,
            pieces = walk.pieces, planned_lengths = planned_lengths)
end

# The length of the planned wire each edge of the path stands for, and the planned length
# of the whole wire. Each piece has its curve length, from the points as given, and its
# staircase, the edges walked between its snapped nodes. A piece whose points snap to one
# node walks no edge: its curve length goes to the next piece that does (around the loop
# for a closed polyline, to the last one at the end of an open polyline). Each walking
# piece's curve length is shared among the edges it walked, in proportion to their
# lengths, so that they add up to it exactly. There is no common factor: edges that cancel
# take their share with them, so the sum falls short of the planned length by what the
# grid cannot hold and nowhere else, and no other part of the wire is stretched for it.
function _planned_lengths(domain::FITDomain, w::PolylineSource, pieces, lengths)
    pts = w.closed ? [w.points..., w.points[1]] : collect(w.points)
    snapped = _polyline_nodes(domain, w)
    at(n) = (domain.nodes_u[n[1]], domain.nodes_v[n[2]], domain.nodes_w[n[3]])
    K = length(pts) - 1
    curve = [hypot((pts[k+1] .- pts[k])...) for k in 1:K]
    stair = [sum(abs, at(snapped[k+1]) .- at(snapped[k])) for k in 1:K]
    walking = findall(>(0), stair)          # not empty: the polyline does not collapse
    covered = zeros(K)                      # the curve length each walking piece stands for
    owner = 0
    for k in K:-1:1                         # each piece to the next walking one
        stair[k] > 0 && (owner = k)
        covered[owner > 0 ? owner : w.closed ? first(walking) : last(walking)] += curve[k]
    end
    return [L * covered[k] / stair[k] for (k, L) in zip(pieces, lengths)], sum(curve)
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

# The walk of a polyline on the grid, kept sparse, with the checks get_source needs:
# - `net`: the current per unit on each edge it passes, ±1, or ±2 where it runs along an
#   edge twice; edges walked as often back as forth are left out;
# - `edges`, `signs`, `nodes`: the path, its steps in order without those back along the
#   edge just walked (see get_source_path), and `pieces`, the piece of the polyline each
#   step was walked for (k: from point k to point k+1).
# Work and memory grow with the length of the path, not with the grid, so that creating a
# polyline, its path and its length stay cheap on large grids; get_source alone needs the
# vector of the whole grid.
function _polyline_walk(config::FITDomain, w::PolylineSource)
    path = _polyline_nodes(config, w)
    all(==(path[1]), path) &&
        throw(ArgumentError("the polyline collapses on the grid: all its points snap to \
                             the same node"))
    node(n) = n[1] + (n[2]-1)*config.Nu + (n[3]-1)*config.Nu*config.Nv
    net = Dict{Int,Float64}()
    edges, signs, nodes, pieces = Int[], Float64[], [node(path[1])], Int[]
    for k in 1:length(path)-1
        _walk(config, path[k], path[k+1]) do edge, sign, reached
            net[edge] = get(net, edge, 0.0) + sign
            if !isempty(edges) && edges[end] == edge && signs[end] == -sign
                # back along the edge just walked
                pop!(edges); pop!(signs); pop!(nodes); pop!(pieces)
            else
                push!(edges, edge); push!(signs, sign); push!(nodes, node(reached)); push!(pieces, k)
            end
        end
    end
    filter!(e -> !iszero(e.second), net)
    isempty(net) &&
        throw(ArgumentError("the polyline cancels itself on the grid: every edge is walked as \
                             often back as forth, so it carries no current"))
    # a closed path can turn back where it closes: its last step undoes its first
    while length(edges) >= 2 && nodes[1] == nodes[end] && edges[1] == edges[end] &&
          signs[1] == -signs[end]
        popfirst!(edges); popfirst!(signs); popfirst!(nodes); popfirst!(pieces)
        pop!(edges); pop!(signs); pop!(nodes); pop!(pieces)
    end

    # Divergence free except at the ends of an open polyline: S̃J = +1 where the current
    # enters, −1 where it leaves, and 0 at every other node the walk touches; an edge
    # takes its current from its lower node to its upper one. True by construction, so,
    # as for the loops, anything else is a bug in the package, not a matter of the input.
    div = _dual_divergence(config, net)
    if path[1] != path[end]
        start, stop = node(path[1]), node(path[end])
        div[start] = get(div, start, 0.0) - 1.0
        div[stop]  = get(div, stop, 0.0) + 1.0
    end
    all(iszero, values(div)) ||
        error("the polyline source is not connected, its divergence is wrong away from the \
               ends: this is a bug in FITToolbox, please report it")
    return (net = net, edges = edges, signs = signs, nodes = nodes, pieces = pieces)
end

# the edges of the walk with their current per unit (see _discretize_source)
_unit_edges(config::FITDomain, w::PolylineSource) = _polyline_walk(config, w).net

# The walk from node P to node Q along the straight line between them: at every step it
# moves one edge along the axis whose next cell-centre plane the line crosses first, so
# it always runs through the nodes nearest to the line. On a tie the lower axis goes
# first. Calls f(edge, sign, node) for each step in order: the edge's index in 1:3Np,
# +1 along its axis or -1 against it, and the node (i, j, k) the step reaches.
function _walk(f, config::FITDomain, P, Q)
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
        cur[axis] += step
        f((axis-1)*Np + ijk[1] + (ijk[2]-1)*Nu + (ijk[3]-1)*Nu*Nv, step, Tuple(cur))
    end
    return nothing
end
