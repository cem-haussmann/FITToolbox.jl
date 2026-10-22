#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 22/10/2026

"""
    get_source_length(domain, name)   -> (planned = …, actual = …)
    get_source_length(domain, obj_id) -> (planned = …, actual = …)

The length of a wire source as planned, i.e. as given to its create function, and as
it is on the grid, in m. The actual length is that of the primal edges its source
vector occupies, each counted as often as the current runs along it:
`Σ |j_e| L_e / |I|` over the edges `e` with length `L_e`, the current `I`.

- polyline (`create_polyline_source!`): `planned` is the sum of the lengths of its
  straight pieces, including the closing one if it is closed; `actual` the length of
  the edge path. Where the path runs back along an edge it came along, e.g. at a
  sharp corner, the two currents cancel and that edge does not count.
- rectangular loop (`create_rectangular_loop_source!`): `planned` and `actual` are
  the two side lengths, `(length_1, length_2)` as given to the create function, i.e.
  along the two in-plane axes in ascending order; on the grid the corners snap to
  nodes.
- circular loop (`create_circular_loop_source!`): `planned` is the circumference
  `2πR`, `actual` the length of the staircase of edges.

A staircase is longer than the line it follows: its length is the sum of the line's
extents along the axes, not its Euclidean length, and refining the grid does not
change that. A circle becomes a staircase of length `8R`, `4/π ≈ 1.27` times its
circumference, on any grid; a diagonal straight piece is longer by up to `√3`. The
staircase follows the line closely all the same, so what it encloses, and with that
its magnetic field, converges.

Throws an `ArgumentError` if no such object exists or it is not one of these sources.

# Example

```julia
create_circular_loop_source!(domain, 0.1, 0.5, 0.5, 0.5, DirZ(); name = "coil")
get_source_length(domain, "coil")      # (planned = 0.628…, actual = 0.8 or close to it)
```
"""
get_source_length(domain::FITDomain, name::AbstractString) =
    get_source_length(domain, _entry(domain, name).obj_id)

function get_source_length(domain::FITDomain, obj_id::Int)
    e = _entry(domain, obj_id)
    src = e.object
    src isa Union{PolylineSource,RectangularLoop,CircularLoop} ||
        throw(ArgumentError("object $obj_id \"$(e.name)\" is a $(nameof(typeof(src))), not a \
                             polyline source, rectangular loop or circular loop"))
    lengths = _edge_lengths_along(domain, _unit_source(domain, src))
    return _source_length(src, lengths)
end

# The length of the edges a unit source occupies, along each axis: Σ |j_e| L_e
function _edge_lengths_along(domain::FITDomain, J)
    Nu, Nv, Np = domain.Nu, domain.Nv, domain.Np
    edges = (domain.edges_u, domain.edges_v, domain.edges_w)
    total = zeros(3)
    for p in findall(!iszero, J)
        axis = (p - 1) ÷ Np + 1
        q = (p - 1) % Np                              # 0-based node index in the block
        ijk = (q % Nu + 1, (q ÷ Nu) % Nv + 1, q ÷ (Nu * Nv) + 1)
        total[axis] += abs(J[p]) * edges[axis][ijk[axis]]
    end
    return total
end

function _source_length(w::PolylineSource, lengths)
    pts = w.closed ? [w.points..., w.points[1]] : collect(w.points)
    planned = sum(hypot((pts[k+1] .- pts[k])...) for k in 1:length(pts)-1)
    return (planned = planned, actual = sum(lengths))
end

# each side twice: once there, once back
function _source_length(r::RectangularLoop, lengths)
    p, q = _inplane_axes(r.normal)
    return (planned = r.lengths, actual = (lengths[p] / 2, lengths[q] / 2))
end

_source_length(c::CircularLoop, lengths) = (planned = 2π * c.radius, actual = sum(lengths))
