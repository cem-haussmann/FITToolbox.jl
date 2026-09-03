#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 11/08/2026
#These functions were generated with Opus 5.0 as an inverse of get_index_entity

function get_position_of_index(config::FITDomain, entity::TopologicalEntity,
                               normal::Normals, i::Integer, j::Integer, k::Integer;
                               units="m")
    (1 <= i <= config.Nu && 1 <= j <= config.Nv && 1 <= k <= config.Nw) ||
        throw(BoundsError("($i, $j, $k) outside 1:$(config.Nu) × 1:$(config.Nv) × 1:$(config.Nw)"))

    x, y, z = _coords_of_ijk(config, entity, normal, i, j, k)
    s = check_units(units)
    return x/s, y/s, z/s
end

function get_position_of_index(config::FITDomain, entity::TopologicalEntity, p::Integer;
                               units="m")
    Nu, Nv, Np = config.Nu, config.Nv, config.Np
    n = _n_components(entity)
    1 <= p <= n*Np || throw(BoundsError("index $p outside 1:$(n*Np) for $entity"))

    normal = (X(), Y(), Z())[(p - 1) ÷ Np + 1]
    q = mod1(p, Np)
    i = mod1(q, Nu)
    j = mod1((q - 1) ÷ Nu + 1, Nv)
    k = (q - 1) ÷ (Nu * Nv) + 1

    return get_position_of_index(config, entity, normal, i, j, k; units)..., normal
end

_n_components(::Union{PrimalNodeDualVolume, DualNodePrimalVolume}) = 1
_n_components(::TopologicalEntity) = 3