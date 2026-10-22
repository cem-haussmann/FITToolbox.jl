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

Throws an `ArgumentError` if `R` or `r` is not positive, if `r ≥ R` (the tube would close
the hole), if the centre lies outside the domain, or if `name` is already used; nothing
is recorded then. Warns if the torus reaches beyond the domain (it is clipped) or fills
no cell.

# Example

```julia
# a ferrite ring core lying flat: centre (50, 50, 30) cm, ring radius 10 cm, tube radius 2 cm
create_torus!(domain, 50, 50, 30, 10, 2, DirZ(); units = "cm", μ_r = 2000.0,
              name = "core", color = "darkgray")
```
"""
function create_torus!(domain::FITDomain, u_c, v_c, w_c, R, r, axis::Direction;
                       units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                       name::String="", color::Union{Nothing,String}=nothing)
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
    σ, ε_r, μ_r = tor.material.σ, tor.material.ε_r, tor.material.μ_r

    R > 0 || throw(ArgumentError("ring radius R must be positive, got $R m"))
    r > 0 || throw(ArgumentError("tube radius r must be positive, got $r m"))
    r < R || throw(ArgumentError("tube radius r = $r m must be smaller than the ring radius \
                                  R = $R m, otherwise the tube closes the hole"))

    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    all(d -> _within(tor.center[d], nodes[d]), 1:3) ||
        throw(ArgumentError("torus centre $(tor.center) m lies outside the domain"))

    a = _axis_index(tor.axis)                 # along the axis
    b, c = filter(!=(a), (1, 2, 3))           # in the ring's plane

    # Reaching past the boundary clips the torus; warn, as for the cylinder
    reach = ntuple(d -> d == a ? r : R + r, 3)
    !all(d -> _within(tor.center[d] - reach[d], nodes[d]) &&
              _within(tor.center[d] + reach[d], nodes[d]), 1:3) &&
        @warn "torus extends beyond the domain and will be clipped" R r center=tor.center

    Nu, Nv, Nw = domain.Nu, domain.Nv, domain.Nw
    md = domain.material

    # Node positions relative to the centre
    du = domain.nodes_u .- tor.center[1]
    dv = domain.nodes_v .- tor.center[2]
    dw = domain.nodes_w .- tor.center[3]
    inside = _inside_test(tor)

    # k outermost for threading (stride Nu*Nv, so threads write far apart),
    # i innermost because md is contiguous in i.
    filled = zeros(Int, Nw-1)
    @threads for k in 1:Nw-1
        ww = (dw[k+1] - dw[k]) / 6.0
        w_samples = (dw[k] + ww, dw[k] + 3.0*ww, dw[k] + 5.0*ww)

        @inbounds for j in 1:Nv-1
            vw = (dv[j+1] - dv[j]) / 6.0
            v_samples = (dv[j] + vw, dv[j] + 3.0*vw, dv[j] + 5.0*vw)

            for i in 1:Nu-1
                uw = (du[i+1] - du[i]) / 6.0
                u_samples = (du[i] + uw, du[i] + 3.0*uw, du[i] + 5.0*uw)

                # 27 sub-voxel centres at 1/6, 3/6, 5/6 of the cell in each direction
                count_inside = 0
                for us in u_samples, vs in v_samples, ws in w_samples
                    inside((us, vs, ws)) && (count_inside += 1)
                end

                # 50% fill rule, as for the sphere: 14 or more of 27 inside
                if count_inside >= 14
                    md[i, j, k, 1] = σ
                    md[i, j, k, 2] = ε_r
                    md[i, j, k, 3] = μ_r
                    filled[k] += 1
                end
            end
        end
    end

    sum(filled) > 0 || @warn "no cells were filled — is the tube radius smaller than one cell?"
end
