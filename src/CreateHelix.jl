#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 22/10/2026

struct Helix <: AbstractSolid
    base::NTuple{3,Float64}    # on the axis, level with the start of the wire
    major_radius::Float64      # R: from the axis to the middle of the wire
    minor_radius::Float64      # r: of the wire, r < R and 2r < pitch
    pitch::Float64             # rise along +axis per turn
    turns::Float64             # may be fractional
    axis::Direction
    start_angle::Float64       # degrees, where the wire starts (see create_helix!)
    left_handed::Bool
    material::Material
end

_describe_geometry(h::Helix) =
    "base=$(_short(h.base)) R=$(_short(h.major_radius)) r=$(_short(h.minor_radius)) " *
    "pitch=$(_short(h.pitch)) turns=$(_short(h.turns)) axis=$(nameof(typeof(h.axis))) " *
    "start=$(_short(h.start_angle))° " * (h.left_handed ? "left-handed" : "right-handed")
_describe_material(h::Helix) = _describe_material(h.material)

"""
    create_helix!(domain, u_b, v_b, w_b, R, r, pitch, turns, axis::Direction;
                  start_angle = 0, left_handed = false, units = "m",
                  σ = 0.0, ε_r = 1.0, μ_r = 1.0, name = "", color = nothing) -> obj_id

Add a helical wire, e.g. the winding of a coil, to `domain` as a solid and return its
`obj_id`. Like `create_cylinder!` and `create_torus!` it only sets the material of the
cells it fills; give it `σ` to make it a conductor.

The wire's centre line winds around `axis` (`DirX()`, `DirY()` or `DirZ()`) at radius
`R`, starting level with `(u_b, v_b, w_b)`, a point on the axis, and rising by `pitch`
along +`axis` per turn for `turns` turns, which may be fractional. The wire has the radius
`r`, measured perpendicular to the wire, and round ends: a half sphere at each end, so
that a lead starting at an end point always connects to it. All lengths are given in
`units`.

`start_angle` (in degrees) is where the wire starts, counted from the in-plane axis that
follows `axis` cyclically (y for `DirX()`, z for `DirY()`, x for `DirZ()`) towards the
other one. The helix is right-handed: it turns counter-clockwise about +`axis` while it
rises, as a loop with that normal; `left_handed = true` turns it the other way. The end
points are returned by `get_helix_ends`.

A cell takes the material when at least half of it lies inside the wire, judged by 27
sample points per cell, the same rule as for the other solids.

With `color = "background"` the helix is a cut-out: it takes the background material of
the domain (`σ`, `ε_r` and `μ_r` are ignored), e.g. to cut a hole into an earlier object,
and `view_domain` shows it as such.

Throws an `ArgumentError` if `R`, `r`, `pitch` or `turns` is not positive, if `r ≥ R`,
if `2r ≥ pitch` (neighbouring turns would touch), if the base lies outside the domain,
or if `name` is already used; nothing is recorded then. Warns if the gap between the
turns, `pitch - 2r`, is smaller than a cell along the axis (the turns can merge on the
grid, which shorts them in a conductor), if the helix reaches beyond the domain (it is
clipped), or if it fills no cell.

# Example

```julia
# a copper winding around a core leg along z: wire radius 6 mm, 5 turns of 20 mm,
# helix radius 35 mm, starting at -x
create_helix!(domain, 0.165, 0.10, 0.15, 0.035, 0.006, 0.02, 5, DirZ();
              start_angle = 180, σ = 5.8e7, name = "coil", color = "copper")
get_helix_ends(domain, "coil")       # where to attach the leads
```
"""
function create_helix!(domain::FITDomain, u_b, v_b, w_b, R, r, pitch, turns, axis::Direction;
                       start_angle::Real=0, left_handed::Bool=false, units="m",
                       σ=0.0, ε_r=1.0, μ_r=1.0,
                       name::AbstractString="", color::Union{Nothing,AbstractString}=nothing)
    unitToMeter = check_units(units)
    helix = Helix((u_b*unitToMeter, v_b*unitToMeter, w_b*unitToMeter),
                  R*unitToMeter, r*unitToMeter, pitch*unitToMeter, turns, axis,
                  start_angle, left_handed, _solid_material(domain, σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, helix)
    return _record_create!(domain, helix; name)
end

"""
    get_helix_ends(domain, name)   -> (start = (u, v, w), stop = (u, v, w))
    get_helix_ends(domain, obj_id) -> (start = (u, v, w), stop = (u, v, w))

The two end points of a helix from `create_helix!`, in m: the centres of its round ends,
on the wire's centre line, where the wire starts and where it stops after its turns.
A lead that starts at one of them is connected to the wire. The ends of a polyline source
come from `get_source_ends` in the same form. Throws an `ArgumentError` if no such object
exists or it is not a helix.
"""
get_helix_ends(domain::FITDomain, name::AbstractString) =
    get_helix_ends(domain, _entry(domain, name).obj_id)

function get_helix_ends(domain::FITDomain, obj_id::Int)
    e = _entry(domain, obj_id)
    e.object isa Helix ||
        throw(ArgumentError("object $obj_id \"$(e.name)\" is a $(nameof(typeof(e.object))), \
                             not a helix"))
    h = e.object
    return (start = _helix_point(h, 0.0), stop = _helix_point(h, 2π * h.turns))
end

# The centre line at parameter t ∈ [0, 2π·turns]: its angle θ0 ± t and height pitch·t/2π
function _helix_point(h::Helix, t)
    a = _axis_index(h.axis)
    b, _ = _cyclic_axes(h.axis)             # the other in-plane axis takes the sine
    θ = deg2rad(h.start_angle) + (h.left_handed ? -t : t)
    off = ntuple(d -> d == a ? h.pitch * t / 2π :
                      d == b ? h.major_radius * cos(θ) : h.major_radius * sin(θ), 3)
    return h.base .+ off
end

# Squared distance from a point to the centre line, given relative to the base in the
# helix's own coordinates: distance ρ from the axis, angle φ (from b towards c) and
# height z. The nearest point lies on one of the two turns closest in height: start at
# the parameter where the line passes the angle φ on each, refine with Newton steps, and
# clamp to the ends, which gives the round caps.
function _helix_distance2(ρ, φ, z, R, P, θ0, s, T)
    f(t) = ρ^2 + R^2 - 2ρ * R * cos(φ - θ0 - s * t) + (z - P * t)^2
    t_angle = mod(s * (φ - θ0), 2π)                  # passes φ at t_angle + 2πk
    k_last = floor((T - t_angle) / 2π)               # the turns that do pass φ: k = 0 … k_last
    k0 = clamp(floor((z / P - t_angle) / 2π), 0, max(k_last - 1, 0))
    # the two turns closest in height that pass φ, and the two ends
    starts = (t_angle + 2π * k0, t_angle + 2π * (k0 + 1), 0.0, T)
    best = Inf
    for t0 in starts
        t = clamp(t0, 0.0, T)
        for _ in 1:4
            ψ = φ - θ0 - s * t
            d1 = -2ρ * R * s * sin(ψ) - 2P * (z - P * t)
            d2 = 2ρ * R * cos(ψ) + 2P^2
            d2 > 0 || break
            t = clamp(t - d1 / d2, 0.0, T)
        end
        best = min(best, f(t))
    end
    return best
end

# The solid as a point test, on coordinates relative to _reference, and its bounding box
_reference(h::Helix) = h.base
function _inside_test(h::Helix)
    a = _axis_index(h.axis)
    b, c = _cyclic_axes(h.axis)
    R, r, height = h.major_radius, h.minor_radius, h.pitch * h.turns
    P, θ0, s, T = h.pitch / 2π, deg2rad(h.start_angle), h.left_handed ? -1.0 : 1.0, 2π * h.turns
    return function (p)
        ρ = hypot(p[b], p[c])
        abs(ρ - R) <= r && -r <= p[a] <= height + r || return false   # the shell around it
        return _helix_distance2(ρ, atan(p[c], p[b]), p[a], R, P, θ0, s, T) <= r^2
    end
end
function _extent(h::Helix)
    a = _axis_index(h.axis)
    R, r = h.major_radius, h.minor_radius
    return (ntuple(d -> d == a ? h.base[d] - r : h.base[d] - R - r, 3),
            ntuple(d -> d == a ? h.base[d] + h.pitch * h.turns + r : h.base[d] + R + r, 3))
end

# Signed distance to the surface, on coordinates relative to _reference: negative inside.
# Only the views use it (smooth surfaces, cut-outs subtracted).
function _signed_distance(h::Helix, p)
    a = _axis_index(h.axis)
    b, c = _cyclic_axes(h.axis)
    d2 = _helix_distance2(hypot(p[b], p[c]), atan(p[c], p[b]), p[a], h.major_radius,
                          h.pitch / 2π, deg2rad(h.start_angle), h.left_handed ? -1.0 : 1.0,
                          2π * h.turns)
    return sqrt(max(d2, 0.0)) - h.minor_radius          # d2 can round to just below 0 on the line
end

function _apply!(domain::FITDomain, hx::Helix)
    R, r, pitch, turns = hx.major_radius, hx.minor_radius, hx.pitch, hx.turns
    σ, ε_r, μ_r = hx.material.σ, hx.material.ε_r, hx.material.μ_r

    R > 0     || throw(ArgumentError("helix radius R must be positive, got $R m"))
    r > 0     || throw(ArgumentError("wire radius r must be positive, got $r m"))
    pitch > 0 || throw(ArgumentError("pitch must be positive, got $pitch m"))
    turns > 0 || throw(ArgumentError("turns must be positive, got $turns"))
    r < R     || throw(ArgumentError("wire radius r = $r m must be smaller than the helix \
                                      radius R = $R m"))
    2r < pitch || throw(ArgumentError("wire radius r = $r m is too thick for the pitch \
                                       $pitch m: neighbouring turns would touch (2r ≥ pitch)"))

    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    all(d -> _within(hx.base[d], nodes[d]), 1:3) ||
        throw(ArgumentError("helix base $(hx.base) m lies outside the domain"))

    a = _axis_index(hx.axis)
    height = pitch * turns

    lo = ntuple(d -> d == a ? hx.base[d] - r : hx.base[d] - R - r, 3)
    hi = ntuple(d -> d == a ? hx.base[d] + height + r : hx.base[d] + R + r, 3)
    !all(d -> _within(lo[d], nodes[d]) && _within(hi[d], nodes[d]), 1:3) &&
        @warn "helix extends beyond the domain and will be clipped" R r pitch turns base=hx.base

    # the turns are apart by pitch - 2r; on a coarser grid they can merge
    span = findall(x -> lo[a] <= x <= hi[a], nodes[a])
    cell = length(span) > 1 ? maximum(diff(nodes[a][first(span):last(span)])) :
                              maximum(diff(nodes[a]))
    pitch - 2r < cell &&
        @warn "the gap between the turns, $(pitch - 2r) m, is smaller than a cell ($cell m): \
               neighbouring turns can merge on the grid" pitch r

    Nu, Nv, Nw = domain.Nu, domain.Nv, domain.Nw
    md = domain.material

    # Node positions relative to the base
    du = domain.nodes_u .- hx.base[1]
    dv = domain.nodes_v .- hx.base[2]
    dw = domain.nodes_w .- hx.base[3]
    inside = _inside_test(hx)

    # k outermost for threading (stride Nu*Nv, so threads write far apart),
    # i innermost because md is contiguous in i.
    # only the cells around the solid can be filled: scan those, not the whole grid
    cells = _cell_ranges(domain, _extent(hx))
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

    sum(filled) > 0 || @warn "no cells were filled — is the wire radius smaller than one cell?"
end
