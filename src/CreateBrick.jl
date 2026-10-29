#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 14/01/2026

struct Brick <: AbstractSolid
    origin::NTuple{3,Float64}
    lengths::NTuple{3,Float64}
    material::Material
end

_describe_geometry(b::Brick) = "origin=$(_short(b.origin)) lengths=$(_short(b.lengths))"
_describe_material(b::Brick) = _describe_material(b.material)

"""
    create_brick!(domain, u_o, v_o, w_o, u_length, v_length, w_length;
                  units = "m", σ = 0.0, ε_r = 1.0, μ_r = 1.0, name = "", color = nothing)
        -> obj_id

Add a solid brick, a box with its sides along the axes, to `domain` and return its
`obj_id`. `(u_o, v_o, w_o)` is its corner with the smallest coordinates; from there it
extends by `u_length`, `v_length` and `w_length`. All lengths are given in `units`. The
object is recorded in the history: `list_objects` shows it, and `remove_object!` and
`undo!` take it away again.

Each face snaps to the nearest node plane, halfway between two to the lower one, as the
loop sources and charges do, and every cell between the faces takes the material: the
brick is exact on the grid wherever its faces fall on nodes.

With `color = "background"` the brick is a cut-out: it takes the background material of
the domain (`σ`, `ε_r` and `μ_r` are ignored), e.g. to cut a window into an earlier
object, and `view_domain` shows it as such.

Throws an `ArgumentError` if a length is not positive, if the corner lies outside the
domain, if the brick has no cell along an axis on the grid (less than about half a cell
of it lies in the domain), or if `name` is already used; nothing is recorded then.
Warns if the brick reaches beyond the domain: like the other solids, it is clipped there.

# Example

```julia
# a copper plate, 80 cm × 80 cm × 2 cm, its corner at (10, 10, 20) cm
create_brick!(domain, 10, 10, 20, 80, 80, 2; units = "cm", σ = 5.8e7, name = "plate",
              color = "copper")
```
"""
function create_brick!(domain::FITDomain, u_o, v_o, w_o, u_length, v_length, w_length;
                       units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                       name::AbstractString="", color::Union{Nothing,AbstractString}=nothing)
    unitToMeter = check_units(units)
    brick = Brick((u_o*unitToMeter, v_o*unitToMeter, w_o*unitToMeter),
                  (u_length*unitToMeter, v_length*unitToMeter, w_length*unitToMeter),
                  _solid_material(domain, σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, brick)
    return _record_create!(domain, brick; name)
end


# The solid as a point test, on coordinates relative to _reference, and its bounding box
# (the material snaps the brick to the nearest nodes instead, see _apply!)
_reference(b::Brick) = b.origin
_inside_test(b::Brick) = (L = b.lengths; p -> 0 <= p[1] <= L[1] && 0 <= p[2] <= L[2] && 0 <= p[3] <= L[3])
_extent(b::Brick) = (b.origin, b.origin .+ b.lengths)
# Signed distance to the surface, on coordinates relative to _reference: negative inside.
# Only the views use it (smooth surfaces, cut-outs subtracted).
function _signed_distance(b::Brick, p)
    q = ntuple(d -> abs(p[d] - b.lengths[d] / 2) - b.lengths[d] / 2, 3)
    return hypot(max(q[1], 0.0), max(q[2], 0.0), max(q[3], 0.0)) + min(max(q[1], q[2], q[3]), 0.0)
end

function _apply!(domain::FITDomain, brick::Brick)
    u_min, v_min, w_min = brick.origin
    u_max, v_max, w_max = brick.origin .+ brick.lengths
    σ, ε_r, μ_r = brick.material.σ, brick.material.ε_r, brick.material.μ_r

    Nodes_U = domain.nodes_u
    Nodes_V = domain.nodes_v
    Nodes_W = domain.nodes_w

    all(>(0), brick.lengths) ||
        throw(ArgumentError("brick lengths must be positive, got $(brick.lengths) m"))
    # The corner must lie inside the grid, as the centre of a sphere must; otherwise the
    # brick is missed entirely or only partly captured, and without this check that goes
    # unnoticed. Check it before snapping: the snapped indices are clamped into range,
    # so testing them cannot detect a brick that lies outside.
    (_within(u_min, Nodes_U) && _within(v_min, Nodes_V) && _within(w_min, Nodes_W)) ||
        throw(ArgumentError("brick corner $(brick.origin) m lies outside the domain"))

    # Each face snaps to the nearest node, as the corners of the loops and the points of
    # the polylines and charges do: halfway between two nodes, to the lower one. A face
    # beyond the last node snaps to it, which clips the brick there.
    i_min, i_max = _find_index(Nodes_U, u_min), _find_index(Nodes_U, u_max)
    j_min, j_max = _find_index(Nodes_V, v_min), _find_index(Nodes_V, v_max)
    k_min, k_max = _find_index(Nodes_W, w_min), _find_index(Nodes_W, w_max)

    # thinner than about half a cell, or clipped to less than that at the boundary
    i_min < i_max ||
        throw(ArgumentError("brick has zero extent in u after snapping to the grid: less than \
                             about half a cell of u = $(_short(u_min)) … $(_short(u_max)) m \
                             lies in the domain"))
    j_min < j_max ||
        throw(ArgumentError("brick has zero extent in v after snapping to the grid: less than \
                             about half a cell of v = $(_short(v_min)) … $(_short(v_max)) m \
                             lies in the domain"))
    k_min < k_max ||
        throw(ArgumentError("brick has zero extent in w after snapping to the grid: less than \
                             about half a cell of w = $(_short(w_min)) … $(_short(w_max)) m \
                             lies in the domain"))

    # Reaching past the boundary clips the brick. That can be intentional — a plate
    # through the whole domain — so warn rather than throw, as for the other solids.
    (_within(u_max, Nodes_U) && _within(v_max, Nodes_V) && _within(w_max, Nodes_W)) ||
        @warn "brick extends beyond the domain and will be clipped" origin=brick.origin lengths=brick.lengths

    Material_distribution = domain.material

    @debug "Creating brick with snapped grid coordinates:" *
    "\n  u: $(Nodes_U[i_min]) m  to  $(Nodes_U[i_max]) m  (requested: $(u_min) m to $(u_max) m)" *
    "\n  v: $(Nodes_V[j_min]) m  to  $(Nodes_V[j_max]) m  (requested: $(v_min) m to $(v_max) m)" *
    "\n  w: $(Nodes_W[k_min]) m  to  $(Nodes_W[k_max]) m  (requested: $(w_min) m to $(w_max) m)"

    #set the properties to the domain
    md = Material_distribution
    @inbounds for i in i_min:i_max-1, j in j_min:j_max-1, k in k_min:k_max-1
        md[i,j,k,1] = σ
        md[i,j,k,2] = ε_r
        md[i,j,k,3] = μ_r
    end
end

@deprecate create_cube!(args...; kwargs...) create_brick!(args...; kwargs...)

