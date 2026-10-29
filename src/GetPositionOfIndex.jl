#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 11/08/2026
#These functions were generated with Opus 5.0 as an inverse of get_index_entity

"""
    get_position_of_index(config, entity, normal, i, j, k; units = "m") -> (x, y, z)
    get_position_of_index(config, entity, p; units = "m")               -> (x, y, z, normal)

The position of an entity of the grid, in `units`: the inverse of `get_index_entity`.
`entity` is `PrimalNode()`, `DualNode()`, `PrimalEdge()`, `DualEdge()`, `PrimalFacet()`,
`DualFacet()`, `PrimalVolume()` or `DualVolume()`. An edge or a facet is at its midpoint,
a node or a volume at its centre; a dual volume is centred on its primal node, a primal
volume (a cell) on its dual node.

The first form takes the grid indices `(i, j, k)` and, for edges and facets, the
component `normal`, `DirX()`, `DirY()` or `DirZ()`; nodes and volumes ignore it. The
second takes `p`, the index in the FIT vector layout, `1:Np` for nodes and volumes and
`1:3Np` for edges and facets, and returns the component as well (`DirX()` for nodes and
volumes).

Throws a `BoundsError` for an index without an entity: outside the grid, or a dead one,
e.g. a u-edge that would leave the last node plane in u.

# Example

```julia
domain = create_domain([1.0, 1.0, 1.0], [0.5, 0.5, 0.5])
p = get_index_entity(domain, PrimalEdge(), DirZ(), 0.5, 0.5, 0.25)    # a w-edge
get_position_of_index(domain, PrimalEdge(), p)                         # (0.5, 0.5, 0.25, DirZ())
```
"""
function get_position_of_index(config::FITDomain, entity::TopologicalEntity,
                               normal::Direction, i::Integer, j::Integer, k::Integer;
                               units="m")
    # N positions along an axis where the entity sits on the nodes, N - 1 where it sits
    # at the edge centres; beyond them there is none
    ranges = map((N, c) -> 1:N - c, (config.Nu, config.Nv, config.Nw), _centred(entity, normal))
    all(map(in, (i, j, k), ranges)) || throw(BoundsError(CartesianIndices(ranges), (i, j, k)))

    x, y, z = _coords_of_ijk(config, entity, normal, i, j, k)
    s = check_units(units)
    return x/s, y/s, z/s
end

function get_position_of_index(config::FITDomain, entity::TopologicalEntity, p::Integer;
                               units="m")
    Nu, Nv, Np = config.Nu, config.Nv, config.Np
    n = _n_components(entity)
    1 <= p <= n*Np || throw(BoundsError(Base.OneTo(n*Np), p))

    normal = (DirX(), DirY(), DirZ())[(p - 1) ÷ Np + 1]
    q = mod1(p, Np)
    i = mod1(q, Nu)
    j = mod1((q - 1) ÷ Nu + 1, Nv)
    k = (q - 1) ÷ (Nu * Nv) + 1

    return get_position_of_index(config, entity, normal, i, j, k; units)..., normal
end

_n_components(::Union{PrimalNodeDualVolume, DualNodePrimalVolume}) = 1
_n_components(::TopologicalEntity) = 3