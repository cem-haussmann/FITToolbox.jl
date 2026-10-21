#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 21/10/2026

struct RectangularLoop{T<:Union{Float64,ComplexF64}} <: AbstractLoopSource
    origin::NTuple{3,Float64}     # the corner with the smallest in-plane coordinates, m
    lengths::NTuple{2,Float64}    # along the two in-plane axes in ascending order, m
    normal::Direction
    current::T                    # A; complex for phasors in the frequency domain
end

# Int/Float32 → Float64, any Complex → ComplexF64
RectangularLoop(origin, lengths, normal, current::Real) =
    RectangularLoop{Float64}(origin, lengths, normal, current)
RectangularLoop(origin, lengths, normal, current::Complex) =
    RectangularLoop{ComplexF64}(origin, lengths, normal, current)

_describe_geometry(r::RectangularLoop) =
    "origin=$(_short(r.origin)) lengths=$(_short(r.lengths)) normal=$(nameof(typeof(r.normal)))"
_describe_material(r::RectangularLoop) = "I=$(_short(r.current)) A"

# The two in-plane axes in ascending order, the order of `lengths`
_inplane_axes(n::Direction) = filter(!=(_axis_index(n)), (1, 2, 3))
# The same axes in cyclic order (b, c), with b × c along +normal
_cyclic_axes(::DirX) = (2, 3)
_cyclic_axes(::DirY) = (3, 1)
_cyclic_axes(::DirZ) = (1, 2)

"""
    create_rectangular_loop_source!(domain, u_o, v_o, w_o, length_1, length_2, normal;
                                    units = "m", current = 1.0, name = "") -> obj_id
    create_rectangular_loop_source!(domain, u_o, v_o, w_o, (length_1, length_2), normal; …)

Add a closed rectangular wire loop carrying `current` (in A) to `domain` and return its
`obj_id`; its discrete current vector comes from `get_source(domain, obj_id)` or
`get_source(domain, name)`, length `3Np` in the edge layout.

As for `create_brick!`, `(u_o, v_o, w_o)` is a corner: the one with the smallest
in-plane coordinates. The loop lies in the plane through it normal to `normal`
(`DirX()`, `DirY()` or `DirZ()`); `length_1` and `length_2` are its sides along the two
in-plane axes in ascending order: v and w for `DirX()`, u and w for `DirY()`, u and v
for `DirZ()`. The lengths can also be given as one tuple. All lengths are in `units`.

The sides run along grid lines, so the loop needs no staircase: the plane and each side
snap to the nearest node, and every edge along the sides carries ±`current`,
counter-clockwise about the normal, i.e. with the magnetic moment along +u, +v or +w for
a positive current. The loop is closed by construction, so its dual divergence
`S̃J = 0` holds exactly, and its moment is exactly `current` times the area of the
snapped rectangle. `current` may be complex; the source vector is then `ComplexF64`.

# Errors

Throws an `ArgumentError` if a side length is not positive, if a side becomes zero on
the grid (shorter than about half a cell), if the plane lies outside
`nodes[1] … nodes[end-2]` along the normal, if the rectangle does not fit within
`nodes[1] … nodes[end-2]` in the plane (the same gap to the boundary as for
`create_circular_loop_source!`), or if `name` is already used. Nothing is recorded then.

# Example

```julia
# 20 cm × 10 cm in the plane z = 30 cm, corner at (40, 45, 30) cm, 1 kA
create_rectangular_loop_source!(domain, 40, 45, 30, 20, 10, DirZ(); units = "cm",
                                current = 1000.0, name = "frame")
J = get_source(domain, "frame")
```
"""
function create_rectangular_loop_source!(domain::FITDomain, u_o, v_o, w_o,
                                         length_1::Real, length_2::Real, normal::Direction;
                                         units="m", current=1.0, name::String="")
    unitToMeter = check_units(units)
    name = _resolve_name(domain, name)
    loop = RectangularLoop((u_o*unitToMeter, v_o*unitToMeter, w_o*unitToMeter),
                           (length_1*unitToMeter, length_2*unitToMeter), normal, current)
    _discretize_source(domain, loop)     # throws on invalid input, before anything is recorded
    return _record_create!(domain, loop; name)
end

function create_rectangular_loop_source!(domain::FITDomain, u_o, v_o, w_o,
                                         lengths::Union{Tuple,AbstractVector},
                                         normal::Direction; kwargs...)
    length(lengths) == 2 ||
        throw(ArgumentError("a rectangle needs 2 side lengths, got $(length(lengths))"))
    return create_rectangular_loop_source!(domain, u_o, v_o, w_o, lengths[1], lengths[2],
                                           normal; kwargs...)
end

function _unit_source(config::FITDomain, r::RectangularLoop)
    nodes = (config.nodes_u, config.nodes_v, config.nodes_w)
    names = ("u", "v", "w")
    a = _axis_index(r.normal)
    p, q = _inplane_axes(r.normal)

    all(>(0), r.lengths) ||
        throw(ArgumentError("side lengths must be positive, got $(r.lengths) m"))
    # The same gap to the boundary as for the circular loop: nodes[1] … nodes[end-2]
    inner(d, x) = _within(x, nodes[d], nodes[d][1], nodes[d][end-2])
    inner(a, r.origin[a]) ||
        throw(ArgumentError("loop plane lies outside the domain in $(names[a])-direction"))
    for (d, L) in ((p, r.lengths[1]), (q, r.lengths[2]))
        inner(d, r.origin[d]) && inner(d, r.origin[d] + L) ||
            throw(ArgumentError("loop does not fit in the domain in $(names[d])-direction"))
    end

    # the plane and the sides on the nearest nodes
    plane = _find_index(nodes[a], r.origin[a])
    lo = Dict(p => _find_index(nodes[p], r.origin[p]), q => _find_index(nodes[q], r.origin[q]))
    hi = Dict(p => _find_index(nodes[p], r.origin[p] + r.lengths[1]),
              q => _find_index(nodes[q], r.origin[q] + r.lengths[2]))
    all(d -> hi[d] > lo[d], (p, q)) ||
        throw(ArgumentError("the loop collapses on the grid: a side of $(r.lengths) m is \
                             shorter than about half a cell"))

    # counter-clockwise about the normal, in the cyclic in-plane axes (b, c):
    # +b along c = c0, +c along b = b1, −b along c = c1, −c along b = b0
    b, c = _cyclic_axes(r.normal)
    Nu, Nv, Np = config.Nu, config.Nv, config.Np
    J = zeros(Float64, 3*Np)
    function edge!(dir, ib, ic, value)
        ijk = ntuple(d -> d == a ? plane : d == b ? ib : ic, 3)
        J[(dir - 1)*Np + ijk[1] + (ijk[2] - 1)*Nu + (ijk[3] - 1)*Nu*Nv] = value
    end
    for m in lo[b]:hi[b]-1
        edge!(b, m, lo[c],  1.0)
        edge!(b, m, hi[c], -1.0)
    end
    for m in lo[c]:hi[c]-1
        edge!(c, hi[b], m,  1.0)
        edge!(c, lo[b], m, -1.0)
    end

    # A current loop must be closed: S̃J = 0, the dual divergence, is discrete current
    # continuity, and the curl–curl system has no solution without it. True by
    # construction; checked as for the circular loop.
    residual = norm(get_divergence(config, Dual()) * J)
    residual > 1e-12 * norm(J) &&
        @warn "The rectangular loop is not divergence free" residual
    return J
end
