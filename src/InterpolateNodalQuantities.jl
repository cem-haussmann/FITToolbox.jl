#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 06/08/2026

function interpolate(domain::FITDomain, ::PrimalNode, datavectorFIT::AbstractVector, x, y, z; units="m")
        x, y, z = convert_to_meter(x, y, z, units)
        result = _findPositionPrimal(domain, x, y, z)
        isnothing(result) && return NaN

        i,j,k = result
        #x-edge is identical because of the offset like nodal values
        val = _GetFieldValue(domain,X(),datavectorFIT,domain.nodes_u ,domain.nodes_v,domain.nodes_w,x, y, z, i, j, k, (i,j,k)->1.0)
        
        return val
end

function interpolate(domain::FITDomain, ::DualNode, datavectorFIT::AbstractVector, x, y, z; units="m")
        x, y, z = convert_to_meter(x, y, z, units)
        result = _findPositionDual(domain, x, y, z)
        isnothing(result) && return NaN

        i_dual,j_dual,k_dual = result
        #x-edge is identical because of the offset like nodal values
        val = _GetFieldValue(domain,X(),datavectorFIT,domain.edges_u_center ,domain.edges_v_center,domain.edges_w_center,x, y, z, i_dual, j_dual, k_dual, (i,j,k)->1.0)
        
        return val
end