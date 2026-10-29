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

Throws an `ArgumentError` if `radius` is not positive, if the centre lies outside the
domain, or if `name` is already used; nothing is recorded then. Warns if the sphere
reaches beyond the domain (it is clipped there) or fills no cell.

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
    σ, ε_r, μ_r = sphere.material.σ, sphere.material.ε_r, sphere.material.μ_r

    r > 0 || throw(ArgumentError("radius must be positive, got $r m"))

    # The centre must lie inside the grid; otherwise the sphere is missed entirely
    # or only partly captured, and without this check that goes unnoticed.
    (_within(start_pos_u, domain.nodes_u) &&
     _within(start_pos_v, domain.nodes_v) &&
     _within(start_pos_w, domain.nodes_w)) ||
        throw(ArgumentError("sphere centre $(sphere.center) m lies outside the domain"))


    # A sphere reaching past the boundary is clipped. That can be intentional —
    # a half sphere at the edge — so warn rather than throw.
    !(_within(start_pos_u - r, domain.nodes_u) && _within(start_pos_u + r, domain.nodes_u) &&
      _within(start_pos_v - r, domain.nodes_v) && _within(start_pos_v + r, domain.nodes_v) &&
      _within(start_pos_w - r, domain.nodes_w) && _within(start_pos_w + r, domain.nodes_w)) &&
        @warn "sphere extends beyond the domain and will be clipped" radius=r centre=(start_pos_u, start_pos_v, start_pos_w)

    Nu, Nv, Nw = domain.Nu, domain.Nv, domain.Nw
    md = domain.material

    # Node positions relative to the sphere centre
    du = domain.nodes_u .- start_pos_u
    dv = domain.nodes_v .- start_pos_v
    dw = domain.nodes_w .- start_pos_w
    inside = _inside_test(sphere)

    # k outermost for threading (stride Nu*Nv, so threads write far apart),
    # i innermost because md is contiguous in i.
    
    # only the cells around the solid can be filled: scan those, not the whole grid
    cells = _cell_ranges(domain, _extent(sphere))
    filled = zeros(Int, Nw-1)
    @threads for k in cells[3]
        ww = (dw[k+1] - dw[k]) / 6.0
        w_samples = (dw[k] + ww, dw[k] + 3.0*ww, dw[k] + 5.0*ww)

        @inbounds for j in cells[2]
            vw = (dv[j+1] - dv[j]) / 6.0
            v_samples = (dv[j] + vw, dv[j] + 3.0*vw, dv[j] + 5.0*vw)

            for i in cells[1]
                uw = (du[i+1] - du[i]) / 6.0
                u_samples = (du[i] + uw, du[i] + 3.0*uw, du[i] + 5.0*uw)

                # 27 sub-voxel centres at 1/6, 3/6, 5/6 of the cell in each direction
                count_inside = 0
                for us in u_samples, vs in v_samples, ws in w_samples
                    inside((us, vs, ws)) && (count_inside += 1)
                end

                # 50% fill rule: assign the material if 14 or more are inside
                if count_inside >= 14
                    md[i, j, k, 1] = σ
                    md[i, j, k, 2] = ε_r
                    md[i, j, k, 3] = μ_r
                    filled[k] += 1
                end
            end
        end
    end

    sum(filled) > 0 || @warn "no cells were filled — is the radius smaller than one cell?"
end