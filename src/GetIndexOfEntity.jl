#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 06/08/2026

using .FieldReconstruction: TopologicalEntity

const PrimalEdgeOrDualFacet = Union{PrimalEdge, DualFacet}
const DualEdgeOrPrimalFacet = Union{DualEdge, PrimalFacet}
const PrimalNodeDualVolume  = Union{PrimalNode, DualVolume}
const DualNodePrimalVolume  = Union{DualNode, PrimalVolume}

const AllNormals = Union{X,Y,Z}

function _find_index(positions::AbstractVector, val::Real)
    idx = clamp(searchsortedlast(positions, val), 1, length(positions))
    
    if idx < length(positions) && abs(positions[idx+1] - val) < abs(positions[idx] - val)
        idx += 1
    end
    
    return idx
end

_get_normal_offset(::X, Np) = 0
_get_normal_offset(::Y, Np) = Np
_get_normal_offset(::Z, Np) = 2 * Np

function get_index_entity(config::FITDomain, entity::PrimalNodeDualVolume, x_pos, y_pos, z_pos; units="m", atol=1e-9)
    return get_index_entity(config, entity, X(), x_pos, y_pos, z_pos; units=units, atol=atol)
end

function get_index_entity(config::FITDomain, entity::DualNodePrimalVolume, x_pos, y_pos, z_pos; units="m", atol=1e-9)
    return get_index_entity(config, entity, X(), x_pos, y_pos, z_pos; units=units, atol=atol)
end

function get_index_entity(config::FITDomain, entity::TopologicalEntity, normals::Normals, x_pos, y_pos, z_pos; units="m", atol=1e-9)
    unitToMeter = check_units(units)
    x = x_pos * unitToMeter
    y = y_pos * unitToMeter
    z = z_pos * unitToMeter

    Np = config.Np

    result = _get_ijk_and_coords(config, entity, normals, x, y, z)
    
    # Check if position was outside the domain
    if isnothing(result)
        return nothing
    end
    
    i, j, k, gx, gy, gz = result

    # Compute Euclidean distance using the exact grid coordinates
    distance = sqrt((x - gx)^2 + (y - gy)^2 + (z - gz)^2)
    
    if distance > atol
        @warn "Topological entity found with distance $(round(distance, sigdigits=4)) m to requested position, larger than set tolerance ($atol)."
    end

    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    offset = _get_normal_offset(normals, Np)
    return 1 + (i-1)*Mu + (j-1)*Mv + (k-1)*Mw + offset
end

####PrimalEdgeOrDualFacet
# X-Normal
function _get_ijk_and_coords(config, ::PrimalEdgeOrDualFacet, ::X, x, y, z)
    i = _find_index(config.edges_u_center, x)
    j = _find_index(config.nodes_v, y)
    k = _find_index(config.nodes_w, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.edges_u_center[i], config.nodes_v[j], config.nodes_w[k]
end

# Y-Normal
function _get_ijk_and_coords(config, ::PrimalEdgeOrDualFacet, ::Y, x, y, z)
    i = _find_index(config.nodes_u, x)
    j = _find_index(config.edges_v_center, y)
    k = _find_index(config.nodes_w, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.nodes_u[i], config.edges_v_center[j], config.nodes_w[k]
end

# Z-Normal
function _get_ijk_and_coords(config, ::PrimalEdgeOrDualFacet, ::Z, x, y, z)
    i = _find_index(config.nodes_u, x)
    j = _find_index(config.nodes_v, y)
    k = _find_index(config.edges_w_center, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.nodes_u[i], config.nodes_v[j], config.edges_w_center[k]
end


####DualEdgeOrPrimalFacet
# X-Normal
function _get_ijk_and_coords(config, ::DualEdgeOrPrimalFacet, ::X, x, y, z)
    i = _find_index(config.nodes_u, x)
    j = _find_index(config.edges_v_center, y)
    k = _find_index(config.edges_w_center, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.nodes_u[i], config.edges_v_center[j], config.edges_w_center[k]
end

# Y-Normal
function _get_ijk_and_coords(config, ::DualEdgeOrPrimalFacet, ::Y, x, y, z)
    i = _find_index(config.edges_u_center, x)
    j = _find_index(config.nodes_v, y)
    k = _find_index(config.edges_w_center, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.edges_u_center[i], config.nodes_v[j], config.edges_w_center[k]
end

# Z-Normal
function _get_ijk_and_coords(config, ::DualEdgeOrPrimalFacet, ::Z, x, y, z)
    i = _find_index(config.edges_u_center, x)
    j = _find_index(config.edges_v_center, y)
    k = _find_index(config.nodes_w, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.edges_u_center[i], config.edges_v_center[j], config.nodes_w[k]
end


####PrimalNodeDualVolume
# All-Normals
function _get_ijk_and_coords(config, ::PrimalNodeDualVolume, ::AllNormals, x, y, z)
    i = _find_index(config.nodes_u, x)
    j = _find_index(config.nodes_v, y)
    k = _find_index(config.nodes_w, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.nodes_u[i], config.nodes_v[j], config.nodes_w[k]
end


####DualNodePrimalVolume
# All-Normals
function _get_ijk_and_coords(config, ::DualNodePrimalVolume, ::AllNormals, x, y, z)
    i = _find_index(config.edges_u_center, x)
    j = _find_index(config.edges_v_center, y)
    k = _find_index(config.edges_w_center, z)
    
    (isnothing(i) || isnothing(j) || isnothing(k)) && return nothing
    
    return i, j, k, config.edges_u_center[i], config.edges_v_center[j], config.edges_w_center[k]
end

# coordinate lookup, used by position_of_index
_coords_of_ijk(c, ::PrimalEdgeOrDualFacet, ::X, i, j, k) =
    (c.edges_u_center[i], c.nodes_v[j],        c.nodes_w[k])
_coords_of_ijk(c, ::PrimalEdgeOrDualFacet, ::Y, i, j, k) =
    (c.nodes_u[i],        c.edges_v_center[j], c.nodes_w[k])
_coords_of_ijk(c, ::PrimalEdgeOrDualFacet, ::Z, i, j, k) =
    (c.nodes_u[i],        c.nodes_v[j],        c.edges_w_center[k])

_coords_of_ijk(c, ::DualEdgeOrPrimalFacet, ::X, i, j, k) =
    (c.nodes_u[i],        c.edges_v_center[j], c.edges_w_center[k])
_coords_of_ijk(c, ::DualEdgeOrPrimalFacet, ::Y, i, j, k) =
    (c.edges_u_center[i], c.nodes_v[j],        c.edges_w_center[k])
_coords_of_ijk(c, ::DualEdgeOrPrimalFacet, ::Z, i, j, k) =
    (c.edges_u_center[i], c.edges_v_center[j], c.nodes_w[k])

_coords_of_ijk(c, ::PrimalNodeDualVolume, ::AllNormals, i, j, k) =
    (c.nodes_u[i],        c.nodes_v[j],        c.nodes_w[k])
_coords_of_ijk(c, ::DualNodePrimalVolume, ::AllNormals, i, j, k) =
    (c.edges_u_center[i], c.edges_v_center[j], c.edges_w_center[k])