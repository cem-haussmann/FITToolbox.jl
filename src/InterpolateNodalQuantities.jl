#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 06/08/2026

function interpolate(domain::FITDomain, ::PrimalNode, datavectorFIT::AbstractVector, x, y, z; units="m")
        _check_length(domain, PrimalNode(), datavectorFIT)
        x, y, z = convert_to_meter(x, y, z, units)
        result = _findPositionPrimal(domain, x, y, z)
        isnothing(result) && return NaN

        i,j,k = result
        #x-edge is identical because of the offset like nodal values
        val = _GetFieldValue(domain,DirX(),datavectorFIT,domain.nodes_u ,domain.nodes_v,domain.nodes_w,x, y, z, i, j, k, (i,j,k)->1.0)
        
        return val
end

function interpolate(domain::FITDomain, ::DualNode, datavectorFIT::AbstractVector, x, y, z; units="m")
        _check_length(domain, DualNode(), datavectorFIT)
        x, y, z = convert_to_meter(x, y, z, units)
        result = _findPositionDual(domain, x, y, z)
        isnothing(result) && return NaN

        i_dual,j_dual,k_dual = result
        #x-edge is identical because of the offset like nodal values
        val = _GetFieldValue(domain,DirX(),datavectorFIT,domain.edges_u_center ,domain.edges_v_center,domain.edges_w_center,x, y, z, i_dual, j_dual, k_dual, (i,j,k)->1.0)

        return val
end

# A volume quantity is integrated over its cells, e.g. a charge: divided by the cell
# volume it is a density, interpolated between the cell centres like nodal values.
# The dual cells are centred on the primal nodes, the primal cells on the dual nodes.
function interpolate(domain::FITDomain, ::DualVolume, datavectorFIT::AbstractVector, x, y, z; units="m")
        _check_length(domain, DualVolume(), datavectorFIT)
        x, y, z = convert_to_meter(x, y, z, units)
        result = _findPositionPrimal(domain, x, y, z)
        isnothing(result) && return NaN

        i,j,k = result
        return _GetFieldValue(domain,DirX(),datavectorFIT,domain.nodes_u,domain.nodes_v,domain.nodes_w,x, y, z, i, j, k,
                              (i,j,k) -> domain.dual_edges_u[i]*domain.dual_edges_v[j]*domain.dual_edges_w[k])
end

function interpolate(domain::FITDomain, ::PrimalVolume, datavectorFIT::AbstractVector, x, y, z; units="m")
        _check_length(domain, PrimalVolume(), datavectorFIT)
        x, y, z = convert_to_meter(x, y, z, units)
        result = _findPositionDual(domain, x, y, z)
        isnothing(result) && return NaN

        i_dual,j_dual,k_dual = result
        return _GetFieldValue(domain,DirX(),datavectorFIT,domain.edges_u_center,domain.edges_v_center,domain.edges_w_center,x, y, z, i_dual, j_dual, k_dual,
                              (i,j,k) -> domain.edges_u[i]*domain.edges_v[j]*domain.edges_w[k])
end