# CreateWire.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: long time ago

struct CircularLoop <: AbstractSource
    center::NTuple{3,Float64}
    radius::Float64
    normal::Direction 
    current::Float64         # A
end

_describe_geometry(c::CircularLoop) =
    "center=$(c.center) radius=$(c.radius) normal=$(nameof(typeof(c.normal)))"
_describe_material(c::CircularLoop) = "I=$(c.current) A"

"""
    create_circular_loop_source(config, radius, center_u, center_v, center_w, normal; units = "m")
        -> Vector{Float64}

Discrete current vector of a closed circular wire loop carrying 1 A, for use as the
right-hand side ĵ of the curl–curl system. The vector has length `3Np` in the edge
layout. Scale it for other currents.

The loop lies in the node plane `u = const` nearest to `center_u`, and only `DirX()`
normals are implemented. On the grid, the circle becomes a staircase: the boundary of
all cells in that plane whose centres lie inside the radius. Every marked v- or w-edge
carries ±1, oriented counter-clockwise about +u, i.e. with the magnetic moment along +u.
The loop is closed by construction, so `GᵀJ = 0` holds exactly. The enclosed area
approaches πR² as the grid is refined.

`radius` and the centre coordinates are given in `units`.

# Errors and warnings

Throws an `ArgumentError` for a `DirY()` or `DirZ()` normal, if `center_u` lies outside
`nodes_u[1] … nodes_u[end-2]`, or if the circle does not fit within
`nodes[1] … nodes[end-2]` in v or w. Warns if the radius is too small to mark any edge.

# Example

```julia
wire = create_circular_loop_source(domain, 50.0, 850.0, 700.0, 700.0, DirX(); units = "mm")
J    = 1000.0 .* wire                   # 1 kA loop
```
"""
function create_circular_loop_source(config::FITDomain, radius, center_u, center_v, center_w,
                                     normal::Direction; units="m", current=1.0, name::String="")
    unitToMeter = check_units(units)
    name = _resolve_name(config, name)
    loop = CircularLoop((center_u*unitToMeter, center_v*unitToMeter, center_w*unitToMeter), radius*unitToMeter, normal, current)
    J = _discretize_source(config, loop) #backup for invalid input, in case there is an error, nothing is recorded
    _record_create!(config, loop; name)
    return J
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