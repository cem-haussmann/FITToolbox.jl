#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 06/10/2026

struct Cylinder <: AbstractSolid
    base::NTuple{3,Float64}    # centre of the base face
    radius::Float64
    height::Float64            # from the base along +axis
    axis::Direction
    material::Material
end

_describe_geometry(c::Cylinder) =
    "base=$(_short(c.base)) radius=$(_short(c.radius)) height=$(_short(c.height)) " *
    "axis=$(nameof(typeof(c.axis)))"
_describe_material(c::Cylinder) = _describe_material(c.material)

_axis_index(::DirX) = 1
_axis_index(::DirY) = 2
_axis_index(::DirZ) = 3

"""
    create_cylinder!(domain, u_b, v_b, w_b, radius, height, axis::Direction;
                     units = "m", σ = 0.0, ε_r = 1.0, μ_r = 1.0, name = "", color = nothing)
        -> obj_id

Add a solid circular cylinder to `domain` and return its `obj_id`. `(u_b, v_b, w_b)` is
the centre of its base face; from there it extends by `height` along +`axis`, which is
`DirX()`, `DirY()` or `DirZ()`. All lengths are given in `units`. The object is recorded
in the history like `create_brick!` and `create_sphere!`.

A cell takes the material when at least half of it lies inside the cylinder, judged by
27 sample points per cell, the same rule as for `create_sphere!`. Along the axis the
cylinder is therefore exact on the grid wherever its end faces fall on nodes; across
it, the circle becomes a staircase that is the same in every cell layer.

With `color = "background"` the cylinder is a cut-out: it takes the background material of
the domain (`σ`, `ε_r` and `μ_r` are ignored), e.g. to cut a hole into an earlier object,
and `view_domain` shows it as such.

Throws an `ArgumentError` if `radius` or `height` is not positive and finite, if the
cylinder lies outside the domain altogether, if the material is not physical (see
`create_domain`), or if `name` is already used; nothing is recorded then. Warns if the
cylinder reaches beyond the domain (it is clipped) or fills no cell. Its base may lie
outside, e.g. for a rod entering through a wall.

# Example

```julia
# a copper rod of 1 cm radius, standing on z = 10 cm, 30 cm long
create_cylinder!(domain, 50, 50, 10, 1, 30, DirZ(); units = "cm", σ = 5.8e7,
                 name = "rod", color = "copper")
```
"""
function create_cylinder!(domain::FITDomain, u_b, v_b, w_b, radius, height, axis::Direction;
                          units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                          name::AbstractString="", color::Union{Nothing,AbstractString}=nothing)
    unitToMeter = check_units(units)
    cylinder = Cylinder((u_b*unitToMeter, v_b*unitToMeter, w_b*unitToMeter),
                        radius*unitToMeter, height*unitToMeter, axis,
                        _solid_material(domain, σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, cylinder)
    return _record_create!(domain, cylinder; name)
end

# The solid as a point test, on coordinates relative to _reference, and its bounding box
_reference(c::Cylinder) = c.base
function _inside_test(c::Cylinder)
    a = _axis_index(c.axis)                   # along the axis
    b, d = filter(!=(a), (1, 2, 3))           # across it
    R, H = c.radius, c.height
    return s -> 0 <= s[a] <= H && hypot(s[b], s[d]) <= R
end
function _extent(c::Cylinder)
    a = _axis_index(c.axis)
    return (ntuple(d -> d == a ? c.base[d] : c.base[d] - c.radius, 3),
            ntuple(d -> d == a ? c.base[d] + c.height : c.base[d] + c.radius, 3))
end

# Signed distance to the surface, on coordinates relative to _reference: negative inside.
# Only the views use it (smooth surfaces, cut-outs subtracted).
function _signed_distance(c::Cylinder, p)
    a = _axis_index(c.axis)
    b, d = filter(!=(a), (1, 2, 3))
    radial = hypot(p[b], p[d]) - c.radius
    axial  = abs(p[a] - c.height / 2) - c.height / 2
    return hypot(max(radial, 0.0), max(axial, 0.0)) + min(max(radial, axial), 0.0)
end

function _apply!(domain::FITDomain, cyl::Cylinder)
    R, H = cyl.radius, cyl.height

    _positive_length(R) || throw(ArgumentError("radius must be positive and finite, got $R m"))
    _positive_length(H) || throw(ArgumentError("height must be positive and finite, got $H m"))

    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    # the base may lie outside, as long as the cylinder reaches into the domain
    _check_reaches_domain(domain, cyl, "cylinder", "base")

    a = _axis_index(cyl.axis)                 # along the axis
    b, c = filter(!=(a), (1, 2, 3))           # across it

    # Reaching past the boundary clips the cylinder. That can be intentional — a rod
    # through the whole domain — so warn rather than throw.
    !(_within(cyl.base[a], nodes[a]) && _within(cyl.base[a] + H, nodes[a]) &&
      all(d -> _within(cyl.base[d] - R, nodes[d]) && _within(cyl.base[d] + R, nodes[d]), (b, c))) &&
        @warn "cylinder extends beyond the domain and will be clipped" radius=R height=H base=cyl.base

    _paint!(domain, cyl) > 0 || @warn "no cells were filled — is the radius or height smaller than one cell?"
end
