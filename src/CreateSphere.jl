#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 14/01/2026

struct Sphere <: AbstractSolid
    center::NTuple{3,Float64}
    radius::Float64
    material::Material
end

_describe_geometry(s::Sphere) = "center=$(_short(s.center)) radius=$(_short(s.radius))"
_describe_material(s::Sphere) = _describe_material(s.material)

"""
    create_sphere!(domain, u_c, v_c, w_c, radius;
                   units = "m", σ = 0.0, ε_r = 1.0, μ_r = 1.0, name = "", color = nothing)
        -> obj_id

Add a solid sphere with its centre at `(u_c, v_c, w_c)` to `domain` and return its
`obj_id`. All lengths are given in `units`. The object is recorded in the history:
`list_objects` shows it, and `remove_object!` and `undo!` take it away again.

A cell takes the material when at least half of it lies inside the sphere, judged by
27 sample points per cell (14 or more inside), so on the grid the sphere becomes a
staircase whose volume approaches `4π/3 radius³` as the grid is refined.
`create_cylinder!`, `create_torus!` and `create_helix!` use the same rule.

With `color = "background"` the sphere is a cut-out: it takes the background material of
the domain (`σ`, `ε_r` and `μ_r` are ignored), e.g. to cut a hole into an earlier
object, and `view_domain` shows it as such.

Throws an `ArgumentError` if `radius` is not positive and finite, if the sphere lies
outside the domain altogether, if the material is not physical (see `create_domain`), or
if `name` is already used; nothing is recorded then. Warns if the sphere reaches beyond
the domain (it is clipped there) or fills no cell. Its centre may lie outside, e.g. for a
cap reaching in through a wall.

# Example

```julia
# a dielectric sphere of 15 cm radius in the middle of a 1 m domain
create_sphere!(domain, 50, 50, 50, 15; units = "cm", ε_r = 4.0, name = "ball")
```
"""
function create_sphere!(domain::FITDomain, u_c, v_c, w_c, radius;
                        units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                        name::AbstractString="", color::Union{Nothing,AbstractString}=nothing)
    unitToMeter = check_units(units)
    sphere = Sphere((u_c*unitToMeter, v_c*unitToMeter, w_c*unitToMeter),
                    radius*unitToMeter,
                    _solid_material(domain, σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, sphere)
    return _record_create!(domain, sphere; name)
end

# The solid as a point test, on coordinates relative to _reference, and its bounding box
_reference(s::Sphere) = s.center
_inside_test(s::Sphere) = (r = s.radius; p -> hypot(p[1], p[2], p[3]) <= r)   # robust against overflow
_extent(s::Sphere) = (s.center .- s.radius, s.center .+ s.radius)
# Signed distance to the surface, on coordinates relative to _reference: negative inside.
# Only the views use it (smooth surfaces, cut-outs subtracted).
_signed_distance(s::Sphere, p) = hypot(p[1], p[2], p[3]) - s.radius

function _apply!(domain::FITDomain, sphere::Sphere)
    start_pos_u, start_pos_v, start_pos_w = sphere.center
    r = sphere.radius

    _positive_length(r) || throw(ArgumentError("radius must be positive and finite, got $r m"))

    # the centre may lie outside, as long as the sphere reaches into the domain
    _check_reaches_domain(domain, sphere, "sphere", "centre")


    # A sphere reaching past the boundary is clipped. That can be intentional —
    # a half sphere at the edge — so warn rather than throw.
    !(_within(start_pos_u - r, domain.nodes_u) && _within(start_pos_u + r, domain.nodes_u) &&
      _within(start_pos_v - r, domain.nodes_v) && _within(start_pos_v + r, domain.nodes_v) &&
      _within(start_pos_w - r, domain.nodes_w) && _within(start_pos_w + r, domain.nodes_w)) &&
        @warn "sphere extends beyond the domain and will be clipped" radius=r centre=(start_pos_u, start_pos_v, start_pos_w)

    _paint!(domain, sphere) > 0 || @warn "no cells were filled — is the radius smaller than one cell?"
end