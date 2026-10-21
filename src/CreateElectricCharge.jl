#This file was created by
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 21/10/2026

struct PointCharge <: AbstractSource
    position::NTuple{3,Float64}   # m, as requested; the grid node is chosen when discretised
    charge::Float64               # C
end

_describe_geometry(q::PointCharge) = "position=$(_short(q.position))"
_describe_material(q::PointCharge) = "q=$(_short(q.charge)) C"

"""
    create_electric_charge!(domain, u, v, w; charge, units = "m", name = "") -> obj_id

Add a point charge at `(u, v, w)` to `domain` and return its `obj_id`. The charge is
recorded in the history like the other objects; its discrete charge vector, the
right-hand side `q` of the Poisson problem `Gᵀ M_ε G Φ = q`, comes from
`get_source(domain, obj_id)` or `get_source(domain, name)`. It has length `Np` and holds
`charge` (in C) in the dual volume of the primal node nearest to the position, zero
elsewhere: charge in FIT is volume-integrated, so the value is the charge itself, not a
density. It is recomputed on the current grid, so it stays valid after
`change_resolution`.

`charge` (in C) is required and may be negative; there is no default, so a charge is
never set without saying so. The position is given in `units`.

Throws an `ArgumentError` if the position lies outside the domain or `name` is already
used; nothing is recorded then.

# Example

```julia
create_electric_charge!(domain, 50.0, 50.0, 50.0; charge = 1.602176634e-19,
                        units = "mm", name = "proton")
q = get_source(domain, "proton")
```
"""
function create_electric_charge!(domain::FITDomain, u, v, w;
                                 charge::Real, units="m", name::String="")
    unitToMeter = check_units(units)
    name = _resolve_name(domain, name)
    pc = PointCharge((u*unitToMeter, v*unitToMeter, w*unitToMeter), charge)
    _discretize_source(domain, pc)     # throws on invalid input, before anything is recorded
    return _record_create!(domain, pc; name)
end

# The primal node nearest to the charge, whose dual volume receives it
function _charge_node(domain::FITDomain, pc::PointCharge)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    all(d -> _within(pc.position[d], nodes[d]), 1:3) ||
        throw(ArgumentError("charge position $(pc.position) m lies outside the domain"))
    return ntuple(d -> _find_index(nodes[d], pc.position[d]), 3)
end

function _discretize_source(domain::FITDomain, pc::PointCharge)
    i, j, k = _charge_node(domain, pc)
    q = zeros(domain.Np)
    q[i + (j-1)*domain.Nu + (k-1)*domain.Nu*domain.Nv] = pc.charge
    return q
end
