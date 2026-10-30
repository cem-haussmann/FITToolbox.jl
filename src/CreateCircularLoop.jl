# CreateCircularLoop.jl (CreateWire.jl until 0.3.0)
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: long time ago

struct CircularLoop{T<:Union{Float64,ComplexF64}} <: AbstractLoopSource
    center::NTuple{3,Float64}
    radius::Float64
    normal::Direction
    current::T               # A; complex for phasors in the frequency domain
end

# Int/Float32 → Float64, any Complex → ComplexF64
CircularLoop(center, radius, normal, current::Real)    = CircularLoop{Float64}(center, radius, normal, current)
CircularLoop(center, radius, normal, current::Complex) = CircularLoop{ComplexF64}(center, radius, normal, current)


_describe_geometry(c::CircularLoop) =
    "center=$(_short(c.center)) radius=$(_short(c.radius)) normal=$(nameof(typeof(c.normal)))"
_describe_material(c::CircularLoop) = "I=$(_short(c.current)) A"

"""
    create_circular_loop_source!(domain, center_u, center_v, center_w, radius, normal;
                                 units = "m", current = 1.0, name = "") -> obj_id

Add a closed circular wire loop carrying `current` (in A) to `domain` and return its
`obj_id`. The loop is recorded in the history like `create_brick!`; its discrete
current vector, the right-hand side ĵ of the curl–curl system, comes from
`get_source(domain, obj_id)` or `get_source(domain, name)`. It has length `3Np` in the
edge layout and is recomputed on the current grid, so it stays valid after
`change_resolution`. `current` may be complex, e.g. for phasors in the frequency
domain; the source vector is then `ComplexF64`.

`normal` is `DirX()`, `DirY()` or `DirZ()`. The loop lies in the node plane normal to
it that is nearest to the centre, e.g. `u = const` nearest to `center_u` for `DirX()`.
On the grid, the circle becomes a staircase: the boundary of all cells in that plane
whose centres lie inside the radius. Every marked in-plane edge carries ±`current`,
oriented counter-clockwise about the normal, i.e. with the magnetic moment along +u,
+v or +w for a positive current. The loop is closed by construction, so its dual
divergence `S̃J = 0` holds exactly. The enclosed area approaches πR² as the grid is
refined.

The centre coordinates and `radius` are given in `units`. The centre comes first, as for
`create_sphere!` and the other objects; the deprecated `create_circular_loop_source`
takes the radius first.

# Errors and warnings

Throws an `ArgumentError` if the centre coordinate along the normal lies outside
`nodes[1] … nodes[end-2]` in that direction, if the circle does not fit within
`nodes[1] … nodes[end-2]` in the two in-plane directions, if the grid has fewer than
three nodes in a direction (that range is then empty), if the radius is not positive
and finite, if `current` is not finite, or if `name` is already used. Nothing is recorded then. Warns if the radius is too
small to mark any edge.

# Example

```julia
create_circular_loop_source!(domain, 850.0, 700.0, 700.0, 50.0, DirX();
                             units = "mm", current = 1000.0, name = "coil")   # 1 kA loop
J = get_source(domain, "coil")
```
"""
function create_circular_loop_source!(domain::FITDomain, center_u, center_v, center_w, radius,
                                     normal::Direction; units="m", current=1.0,
                                     name::AbstractString="")
    unitToMeter = check_units(units)
    name = _resolve_name(domain, name)
    loop = CircularLoop((center_u*unitToMeter, center_v*unitToMeter, center_w*unitToMeter), radius*unitToMeter, normal, current)
    _check_source(domain, loop)          # throws on invalid input, before anything is recorded
    return _record_create!(domain, loop; name)
end

# Loops lie within nodes[1] … nodes[end-2] in every direction (see the docstring), so the
# grid needs at least three nodes in each; the rectangular loop shares this
function _check_loop_grid(config::FITDomain)
    for (name, n) in (("u", config.Nu), ("v", config.Nv), ("w", config.Nw))
        n >= 3 || throw(ArgumentError("a loop source needs at least 3 nodes in every \
                                       direction, but the grid has $n in $name"))
    end
    return nothing
end

# The checks of the docstring: a positive radius, the plane within nodes[1] … nodes[end-2]
# along the normal, and the circle within them in the two in-plane directions
function _check_loop_fits(config::FITDomain, loop::CircularLoop)
    R = loop.radius
    _positive_length(R) || throw(ArgumentError("radius must be positive and finite, got $R m"))
    _check_loop_grid(config)
    nodes = (config.nodes_u, config.nodes_v, config.nodes_w)
    names = ("u", "v", "w")
    inner(d, x) = _within(x, nodes[d], nodes[d][1], nodes[d][end-2])
    a = _axis_index(loop.normal)
    inner(a, loop.center[a]) ||
        throw(ArgumentError("the coil plane does not fit in the domain in $(names[a])-direction: \
                             it must lie within $(_loop_range(nodes, names, a)), got \
                             $(names[a]) = $(_short(loop.center[a])) m"))
    for d in _inplane_axes(loop.normal)
        inner(d, loop.center[d] - R) && inner(d, loop.center[d] + R) ||
            throw(ArgumentError("coil does not fit in the domain in $(names[d])-direction: \
                                 it must lie within $(_loop_range(nodes, names, d)), got \
                                 $(_short(loop.center[d] - R)) … $(_short(loop.center[d] + R)) m"))
    end
    return nothing
end

# The range a loop must lie within along axis d, for the error messages of both loops
_loop_range(nodes, names, d) =
    "$(names[d]) = $(_short(nodes[d][1])) … $(_short(nodes[d][end-2])) m (nodes[1] … nodes[end-2])"

# A current loop must be closed: S̃J = 0, the dual divergence, is discrete current
# continuity, and the curl–curl system has no solution without it. The loops are closed
# by construction, so a residual is a bug in the package, not a matter of the input.
# `edges` holds the loop's current per unit, edge => ±1, so the sums are exact.
function _check_closed(config::FITDomain, edges, what)
    all(iszero, values(_dual_divergence(config, edges))) ||
        error("the $what is not divergence free: this is a bug in FITToolbox, please \
               report it")
    return nothing
end

# The marked edges are the boundary of the cells in the loop's plane whose centres lie
# inside the radius. An edge on the first node line of a direction has no cell before
# it; that cell counts as outside, so a loop reaching into the first cells is closed
# there too. The last cells never lie inside: the circle ends at nodes[end-2].
# Only the cells around the circle are visited, so work and memory grow with its size
# on the grid, not with the plane.
function _unit_edges(config::FITDomain, loop::CircularLoop)
    _check_loop_fits(config, loop)
    nodes   = (config.nodes_u, config.nodes_v, config.nodes_w)
    centres = (config.edges_u_center, config.edges_v_center, config.edges_w_center)
    Nu, Nv, Np = config.Nu, config.Nv, config.Np
    a = _axis_index(loop.normal)
    b, c = _cyclic_axes(loop.normal)                # b × c along +normal
    x, R = loop.center, loop.radius
    plane = _find_index(nodes[a], x[a])

    # The cells along axis d whose centres lie within R of the centre, and two more on
    # each side: the cells just outside decide the edges around the circle, the next ones
    # guard against rounding. Cells outside the circle mark no edge, so more do no harm.
    function cells(d)
        n = length(centres[d])
        return max(1, searchsortedlast(centres[d], x[d] - R) - 1):min(n, searchsortedfirst(centres[d], x[d] + R) + 1)
    end
    # cell (ib, ic) in the plane, by its indices along b and c; cell 0 is before the first
    outside(ib, ic) = ib == 0 || ic == 0 ||
                      hypot(centres[b][ib] - x[b], centres[c][ic] - x[c]) >= R
    # the edge along axis `dir` starting from node (ib, ic) in the plane
    function edge(dir, ib, ic)
        i, j, k = ntuple(d -> d == a ? plane : d == b ? ib : ic, 3)
        return (dir - 1)*Np + i + (j - 1)*Nu + (k - 1)*Nu*Nv
    end

    net = Dict{Int,Float64}()
    for ic in cells(c), ib in cells(b)
        # The edges from node (ib, ic) separate cell (ib, ic) from the cell before it
        # along c (the b-edge) or along b (the c-edge); an edge is marked when exactly
        # one of the two lies inside. Counter-clockwise about the normal, i.e. along
        # normal × r, the direction follows from which side of the edge is inside (the
        # side of the loop it is on), so the loop closes on any grid, also where it
        # collapses onto the centre line:
        inside = !outside(ib, ic)
        # b-edge: +b at the bottom of the loop (inside after it in c), −b at the top
        outside(ib, ic - 1) == inside && (net[edge(b, ib, ic)] = inside ? 1.0 : -1.0)
        # c-edge: −c on the near side in b (inside after it), +c on the far side
        outside(ib - 1, ic) == inside && (net[edge(c, ib, ic)] = inside ? -1.0 : 1.0)
    end

    isempty(net) &&
        @warn "No edges were marked for the wire — is the radius smaller than one cell?"
    _check_closed(config, net, "circular coil")
    return net
end

# deprecated in 0.3.0: the loop is recorded in the history and returns its obj_id, and
# the centre comes before the radius, as for the other objects. The old name keeps the
# old order and still returns the source vector.
@deprecate(create_circular_loop_source(domain::FITDomain, radius, center_u, center_v, center_w,
                                       normal; kwargs...),
           get_source(domain, create_circular_loop_source!(domain, center_u, center_v, center_w,
                                                           radius, normal; kwargs...)))
