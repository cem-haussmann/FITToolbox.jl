#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 22/10/2026

struct Torus <: AbstractSolid
    center::NTuple{3,Float64}  # centre of the ring, the middle of the hole
    major_radius::Float64      # R: from the centre to the middle of the tube
    minor_radius::Float64      # r: of the tube, r < R
    axis::Direction            # symmetry axis, the normal of the ring's plane
    material::Material
end

_describe_geometry(t::Torus) =
    "center=$(_short(t.center)) R=$(_short(t.major_radius)) r=$(_short(t.minor_radius)) " *
    "axis=$(nameof(typeof(t.axis)))"
_describe_material(t::Torus) = _describe_material(t.material)

"""
    create_torus!(domain, u_c, v_c, w_c, R, r, axis::Direction;
                  units = "m", σ = 0.0, ε_r = 1.0, μ_r = 1.0, name = "", color = nothing)
        -> obj_id

Add a solid torus (a ring) to `domain` and return its `obj_id`. `(u_c, v_c, w_c)` is
the centre of the ring, the middle of its hole; `R` is the ring radius, from the centre
to the middle of the tube, and `r` the radius of the tube. `axis`, `DirX()`, `DirY()` or
`DirZ()`, is its symmetry axis, the normal of the ring's plane, as for
`create_circular_loop_source!`: `DirZ()` lays the ring flat in the x–y plane. All lengths
are given in `units`. The object is recorded in the history like the other solids.

A cell takes the material when at least half of it lies inside the torus, judged by 27
sample points per cell, the same rule as for `create_sphere!` and `create_cylinder!`.

With `color = "background"` the torus is a cut-out: it takes the background material of
the domain (`σ`, `ε_r` and `μ_r` are ignored), e.g. to cut a hole into an earlier object,
and `view_domain` shows it as such.

Throws an `ArgumentError` if `R` or `r` is not positive and finite, if `r ≥ R` (the tube
would close the hole), if the torus lies outside the domain altogether, if the material
is not physical (see `create_domain`), or if `name` is already used; nothing is recorded
then. Warns if the torus reaches beyond the domain (it is clipped) or fills no cell. Its
centre may lie outside, as long as part of the ring reaches in.

# Example

```julia
# a ferrite ring core lying flat: centre (50, 50, 30) cm, ring radius 10 cm, tube radius 2 cm
create_torus!(domain, 50, 50, 30, 10, 2, DirZ(); units = "cm", μ_r = 2000.0,
              name = "core", color = "darkgray")
```
"""
function create_torus!(domain::FITDomain, u_c, v_c, w_c, R, r, axis::Direction;
                       units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                       name::AbstractString="", color::Union{Nothing,AbstractString}=nothing)
    unitToMeter = check_units(units)
    torus = Torus((u_c*unitToMeter, v_c*unitToMeter, w_c*unitToMeter),
                  R*unitToMeter, r*unitToMeter, axis, _solid_material(domain, σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, torus)
    return _record_create!(domain, torus; name)
end

# The solid as a point test, on coordinates relative to _reference, and its bounding box
_reference(t::Torus) = t.center
function _inside_test(t::Torus)
    a = _axis_index(t.axis)                   # along the axis
    b, c = filter(!=(a), (1, 2, 3))           # in the ring's plane
    R, r = t.major_radius, t.minor_radius
    # the distance from the tube's centre circle is at most r
    return s -> (hypot(s[b], s[c]) - R)^2 + s[a]^2 <= r^2
end
function _extent(t::Torus)
    a = _axis_index(t.axis)
    reach = ntuple(d -> d == a ? t.minor_radius : t.major_radius + t.minor_radius, 3)
    return (t.center .- reach, t.center .+ reach)
end

# Signed distance to the surface, on coordinates relative to _reference: negative inside.
# Only the views use it (smooth surfaces, cut-outs subtracted).
function _signed_distance(t::Torus, p)
    a = _axis_index(t.axis)
    b, c = filter(!=(a), (1, 2, 3))
    return hypot(hypot(p[b], p[c]) - t.major_radius, p[a]) - t.minor_radius
end

function _apply!(domain::FITDomain, tor::Torus)
    R, r = tor.major_radius, tor.minor_radius

    _positive_length(R) || throw(ArgumentError("ring radius R must be positive and finite, got $R m"))
    _positive_length(r) || throw(ArgumentError("tube radius r must be positive and finite, got $r m"))
    r < R || throw(ArgumentError("tube radius r = $r m must be smaller than the ring radius \
                                  R = $R m, otherwise the tube closes the hole"))

    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    # the centre may lie outside, as long as the torus reaches into the domain
    _check_reaches_domain(domain, tor, "torus", "centre")

    a = _axis_index(tor.axis)                 # along the axis

    # Reaching past the boundary clips the torus; warn, as for the cylinder
    reach = ntuple(d -> d == a ? r : R + r, 3)
    !all(d -> _within(tor.center[d] - reach[d], nodes[d]) &&
              _within(tor.center[d] + reach[d], nodes[d]), 1:3) &&
        @warn "torus extends beyond the domain and will be clipped" R r center=tor.center

    _paint!(domain, tor) > 0 || @warn "no cells were filled — is the tube radius smaller than one cell?"
end
