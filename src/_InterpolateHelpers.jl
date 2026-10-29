#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 21/03/2026

# One value per entity: Np for the nodal and volume quantities, 3Np for edges and facets.
# Without this check a vector of the other kind is read silently, e.g. the u-components
# of an edge vector as nodal values.
_n_values(domain::FITDomain, ::Union{PrimalNode, DualNode, PrimalVolume, DualVolume}) = domain.Np
_n_values(domain::FITDomain, ::Union{PrimalEdge, DualEdge, PrimalFacet, DualFacet}) = 3 * domain.Np

function _check_length(domain::FITDomain, entity, data::AbstractVector)
    length(data) == _n_values(domain, entity) ||
        throw(DimensionMismatch("$(nameof(typeof(entity))) data has length $(length(data)), \
                                 the domain needs $(_n_values(domain, entity))"))
    return nothing
end

function _find_dual_index(centers::AbstractVector, val::Real)
    n = length(centers)
    # A single cell along this axis has a single dual node: the value is constant along
    # it (_GetFieldValue has no neighbour to interpolate to)
    n == 1 && return 1

    # Within half a cell of the boundary the point lies beyond the first or last dual
    # node: the outermost pair is used and the value extrapolated over at most half a
    # cell. Points outside the domain never get here (see _findPositionDual).
    return clamp(searchsortedlast(centers, val), 1, n - 1)
end

function _find_primal_index(nodes::AbstractVector, val::Real)
    if !_within(val, nodes)
        @warn "The defined position ($val) is outside of the domain!"
        return nothing
    end
    return clamp(searchsortedlast(nodes, val), 1, length(nodes) - 1)
end

function _findPositionPrimal(domain::FITDomain, x, y, z)
    i = _find_primal_index(domain.nodes_u, x)
    j = _find_primal_index(domain.nodes_v, y)
    k = _find_primal_index(domain.nodes_w, z)

    # If any coordinate is outside the domain, return nothing
    if isnothing(i) || isnothing(j) || isnothing(k)
        return nothing
    end

    return i, j, k
end

function _findPositionDual(domain::FITDomain, x, y, z)
    # Outside the domain there is no value, as for the primal quantities: nothing,
    # which the callers return as NaN (_findPositionPrimal warns).
    isnothing(_findPositionPrimal(domain, x, y, z)) && return nothing
    # Inside it, a point in the outermost half cell has a dual node on one side only;
    # _find_dual_index then takes the outermost pair and the value is extrapolated.
    i_dual = _find_dual_index(domain.edges_u_center, x)
    j_dual = _find_dual_index(domain.edges_v_center, y)
    k_dual = _find_dual_index(domain.edges_w_center, z)

    return i_dual, j_dual, k_dual
end

function _GetFieldValue(domain::FITDomain,normals::Direction,edgeQuantities::AbstractVector,x_pos_array,y_pos_array,z_pos_array,x,y,z,i,j,k,normalization_)
    offset = 0
    if normals isa DirY
        offset+=domain.Np
    elseif normals isa DirZ
        offset+=2*domain.Np
    end

    Nu = domain.Nu
    Nv = domain.Nv

    Np = domain.Np
    
    Mu = 1
    Mv = Nu
    Mw = Nu*Nv
    
    # The step to the next position along each axis: 1, or 0 along an axis with a single
    # position (the one dual node of a single cell), where the value is constant
    di = i < length(x_pos_array) ? 1 : 0
    dj = j < length(y_pos_array) ? 1 : 0
    dk = k < length(z_pos_array) ? 1 : 0

    x_d = di == 0 ? 0.0 : (x-x_pos_array[i])/(x_pos_array[i+1]-x_pos_array[i])
    y_d = dj == 0 ? 0.0 : (y-y_pos_array[j])/(y_pos_array[j+1]-y_pos_array[j])
    z_d = dk == 0 ? 0.0 : (z-z_pos_array[k])/(z_pos_array[k+1]-z_pos_array[k])

    p0 = 1 + (i - 1) * Mu + (j - 1) * Mv + (k - 1) * Mw + offset;
    su, sv, sw = di * Mu, dj * Mv, dk * Mw

    # +1 in u-direction
    p1 = p0 + su;
    # +1 in v-direction
    p2 = p0 + sv;
    # +1 in v-direction; +1 in u direction
    p3 = p0 + sv + su;
    # +1 in w-direction
    p4 = p0 + sw;
    # +1 in w-direction; +1 in u-direction
    p5 = p0 + sw + su;
    # +1 in w-direction; +1 in v-direction
    p6 = p0 + sw + sv;
    # +1 in w-direction; +1 in v-direction; +1 in u-direction
    p7 = p0 + sw + sv + su;

    c00 = edgeQuantities[p0]/(normalization_(i,j,k))*(1-x_d) + edgeQuantities[p1]/(normalization_(i+di,j,k))*x_d
    c01 = edgeQuantities[p4]/(normalization_(i,j,k+dk))*(1-x_d) + edgeQuantities[p5]/(normalization_(i+di,j,k+dk))*x_d
    c10 = edgeQuantities[p2]/(normalization_(i,j+dj,k))*(1-x_d) + edgeQuantities[p3]/(normalization_(i+di,j+dj,k))*x_d
    c11 = edgeQuantities[p6]/(normalization_(i,j+dj,k+dk))*(1-x_d) + edgeQuantities[p7]/(normalization_(i+di,j+dj,k+dk))*x_d

    c0 = c00*(1-y_d) + c10*y_d
    c1 = c01*(1-y_d) + c11*y_d

    c = c0*(1-z_d) + c1*z_d    
    return c
end