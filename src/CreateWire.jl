# CreateWire.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: long time ago

struct CircularLoop{T<:Union{Float64,ComplexF64}} <: AbstractSource
    center::NTuple{3,Float64}
    radius::Float64
    normal::Direction
    current::T               # A; complex for phasors in the frequency domain
end

# Int/Float32 → Float64, any Complex → ComplexF64
CircularLoop(center, radius, normal, current::Real)    = CircularLoop{Float64}(center, radius, normal, current)
CircularLoop(center, radius, normal, current::Complex) = CircularLoop{ComplexF64}(center, radius, normal, current)


_describe_geometry(c::CircularLoop) =
    "center=$(_short(c.center)) radius=$(_short(c.radius)) normal=$(nameof(typeof(c.normal)))"
_describe_material(c::CircularLoop) = "I=$(_short(c.current)) A"

"""
    create_circular_loop_source!(config, radius, center_u, center_v, center_w, normal;
                                 units = "m", current = 1.0, name = "") -> obj_id

Add a closed circular wire loop carrying `current` (in A) to `config` and return its
`obj_id`. The loop is recorded in the history like `create_brick!`; its discrete
current vector, the right-hand side ĵ of the curl–curl system, comes from
`get_source(config, obj_id)` or `get_source(config, name)`. It has length `3Np` in the
edge layout and is recomputed on the current grid, so it stays valid after
`change_resolution`. `current` may be complex, e.g. for phasors in the frequency
domain; the source vector is then `ComplexF64`.

The loop lies in the node plane `u = const` nearest to `center_u`, and only `DirX()`
normals are implemented. On the grid, the circle becomes a staircase: the boundary of
all cells in that plane whose centres lie inside the radius. Every marked v- or w-edge
carries ±`current`, oriented counter-clockwise about +u, i.e. with the magnetic moment
along +u for a positive current. The loop is closed by construction, so `GᵀJ = 0` holds
exactly. The enclosed area approaches πR² as the grid is refined.

`radius` and the centre coordinates are given in `units`.

# Errors and warnings

Throws an `ArgumentError` for a `DirY()` or `DirZ()` normal, if `center_u` lies outside
`nodes_u[1] … nodes_u[end-2]`, if the circle does not fit within
`nodes[1] … nodes[end-2]` in v or w, or if `name` is already used. Nothing is recorded
then. Warns if the radius is too small to mark any edge.

# Example

```julia
create_circular_loop_source!(domain, 50.0, 850.0, 700.0, 700.0, DirX();
                             units = "mm", current = 1000.0, name = "coil")   # 1 kA loop
J = get_source(domain, "coil")
```
"""
function create_circular_loop_source!(config::FITDomain, radius, center_u, center_v, center_w,
                                     normal::Direction; units="m", current=1.0, name::String="")
    unitToMeter = check_units(units)
    name = _resolve_name(config, name)
    loop = CircularLoop((center_u*unitToMeter, center_v*unitToMeter, center_w*unitToMeter), radius*unitToMeter, normal, current)
    _discretize_source(config, loop)     # throws on invalid input, before anything is recorded
    return _record_create!(config, loop; name)
end

_discretize_source(config::FITDomain, loop::CircularLoop) =
    loop.current .* _circular_loop_unit_source(config, loop, loop.normal)

_circular_loop_unit_source(config, loop::CircularLoop, ::DirY) =
    throw(ArgumentError("only DirX() normals are implemented; got DirY()"))

_circular_loop_unit_source(config, loop::CircularLoop, ::DirZ) =
    throw(ArgumentError("only DirX() normals are implemented; got DirZ()"))

function _circular_loop_unit_source(config, loop::CircularLoop, ::DirX)
    Nodes_V, Nodes_W = config.nodes_v, config.nodes_w
    Ev_c, Ew_c = config.edges_v_center, config.edges_w_center
    Nu, Nv, Nw, Np = config.Nu, config.Nv, config.Nw, config.Np

    cu, cv, cw = loop.center
    R = loop.radius
    R > 0 || throw(ArgumentError("radius must be positive, got $R m"))
    
    # The loop lies in the v–w plane, so only v and w constrain the radius; u only
    # has to fall inside the domain.
    config.nodes_u[1] <= cu <= config.nodes_u[end-2] ||
        throw(ArgumentError("coil plane lies outside the domain in u-direction"))
    cv - R >= Nodes_V[1] && cv + R <= Nodes_V[end-2] ||
        throw(ArgumentError("coil does not fit in the domain in v-direction"))
    cw - R >= Nodes_W[1] && cw + R <= Nodes_W[end-2] ||
        throw(ArgumentError("coil does not fit in the domain in w-direction"))
 

    i = _find_index(config.nodes_u, cu)
    
    J = zeros(Float64, 3*Np)
    ρ(v, w) = hypot(v - cv, w - cw)          # distance from the loop centre
 
    for kI in 2:Nw-1, jI in 2:Nv-1
        p = 1 + (i-1) + (jI-1)*Nu + (kI-1)*Nu*Nv
 
        # Three corners of the cell in the v–w plane. An edge is crossed when the
        # circle separates its two endpoints, i.e. exactly one of them is outside.
        out1 = ρ(Ev_c[jI],   Ew_c[kI-1]) >= R
        out2 = ρ(Ev_c[jI-1], Ew_c[kI]  ) >= R
        out3 = ρ(Ev_c[jI],   Ew_c[kI]  ) >= R
 
        # v-directed edge: sign flips across the loop centre in w
        out3 != out1 && (J[p + Np]   = Nodes_W[kI] >= cw ? -1.0 : 1.0)
        # w-directed edge: sign flips across the loop centre in v
        out3 != out2 && (J[p + 2*Np] = Nodes_V[jI] >= cv ?  1.0 : -1.0)
    end
 
    any(!iszero, J) ||
        @warn "No edges were marked for the wire — is the radius smaller than one cell?"
 
    # A current loop must be closed: GᵀJ = 0 is discrete current continuity, and
    # the curl–curl system has no solution without it.
    residual = norm(transpose(get_gradient(config, Primal())) * J)
    residual > 1e-12 * norm(J) &&
        @warn "The circular coil is not divergence free" residual
 
    return J
end

# deprecated in 0.3.0: the loop is recorded in the history and returns its obj_id
function create_circular_loop_source(config::FITDomain, args...; kwargs...)
    Base.depwarn("`create_circular_loop_source` is deprecated, use \
                  `id = create_circular_loop_source!(...)` and `J = get_source(domain, id)`",
                 :create_circular_loop_source)
    return get_source(config, create_circular_loop_source!(config, args...; kwargs...))
end
