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

Throws an `ArgumentError` if `radius` or `height` is not positive, if the centre of the
base face lies outside the domain, or if `name` is already used; nothing is recorded
then. Warns if the cylinder reaches beyond the domain (it is clipped) or fills no cell.

# Example

```julia
# a copper rod of 1 cm radius, standing on z = 10 cm, 30 cm long
create_cylinder!(domain, 50, 50, 10, 1, 30, DirZ(); units = "cm", σ = 5.8e7,
                 name = "rod", color = "copper")
```
"""
function create_cylinder!(domain::FITDomain, u_b, v_b, w_b, radius, height, axis::Direction;
                          units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                          name::String="", color::Union{Nothing,String}=nothing)
    unitToMeter = check_units(units)
    cylinder = Cylinder((u_b*unitToMeter, v_b*unitToMeter, w_b*unitToMeter),
                        radius*unitToMeter, height*unitToMeter, axis,
                        Material(σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, cylinder)
    return _record_create!(domain, cylinder; name)
end

function _apply!(domain::FITDomain, cyl::Cylinder)
    R, H = cyl.radius, cyl.height
    σ, ε_r, μ_r = cyl.material.σ, cyl.material.ε_r, cyl.material.μ_r

    R > 0 || throw(ArgumentError("radius must be positive, got $R m"))
    H > 0 || throw(ArgumentError("height must be positive, got $H m"))

    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    # The base must lie inside the grid; otherwise the cylinder is missed entirely
    # or only partly captured, and without this check that goes unnoticed.
    all(d -> nodes[d][1] <= cyl.base[d] <= nodes[d][end], 1:3) ||
        throw(ArgumentError("cylinder base $(cyl.base) m lies outside the domain"))

    a = _axis_index(cyl.axis)                 # along the axis
    b, c = filter(!=(a), (1, 2, 3))           # across it

    # Reaching past the boundary clips the cylinder. That can be intentional — a rod
    # through the whole domain — so warn rather than throw.
    (cyl.base[a] + H > nodes[a][end] ||
     any(d -> cyl.base[d] - R < nodes[d][1] || cyl.base[d] + R > nodes[d][end], (b, c))) &&
        @warn "cylinder extends beyond the domain and will be clipped" radius=R height=H base=cyl.base

    Nu, Nv, Nw = domain.Nu, domain.Nv, domain.Nw
    md = domain.material

    # Node positions relative to the centre of the base face
    du = domain.nodes_u .- cyl.base[1]
    dv = domain.nodes_v .- cyl.base[2]
    dw = domain.nodes_w .- cyl.base[3]
    inside(s) = 0 <= s[a] <= H && hypot(s[b], s[c]) <= R

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

    sum(filled) > 0 || @warn "no cells were filled — is the radius or height smaller than one cell?"
end
