# InterpolateDualEdgePrimalFacet.jl
function interpolate(domain::FITDomain, ::DualEdge, datavectorFIT::AbstractVector, x, y, z; units="m")
        x, y, z = convert_to_meter(x, y, z, units)

        result = _findPositionPrimal(domain, x, y, z)
        isnothing(result) && return NaN, NaN, NaN

        i,j,k = result
        result = _findPositionDual(domain, x, y, z)
        isnothing(result) && return NaN, NaN, NaN

        i_dual,j_dual,k_dual = result

        Ex = _GetFieldValue(domain,X(), datavectorFIT,domain.nodes_u ,domain.edges_v_center,domain.edges_w_center,x, y, z, i, j_dual, k_dual, (i, j, k) -> domain.dual_edges_u[i])
        Ey = _GetFieldValue(domain,Y(), datavectorFIT,domain.edges_u_center ,domain.nodes_v,domain.edges_w_center,x, y, z, i_dual, j, k_dual, (i, j, k) -> domain.dual_edges_v[j])
        Ez = _GetFieldValue(domain,Z(), datavectorFIT,domain.edges_u_center ,domain.edges_v_center,domain.nodes_w,x, y, z, i_dual, j_dual, k, (i, j, k) -> domain.dual_edges_w[k])
        return Ex, Ey, Ez
end

function interpolate(domain::FITDomain, ::PrimalFacet, datavectorFIT::AbstractVector, x, y, z; units="m")
        x, y, z = convert_to_meter(x, y, z, units)

        result = _findPositionPrimal(domain, x, y, z)
        isnothing(result) && return NaN, NaN, NaN

        i,j,k = result
        result = _findPositionDual(domain, x, y, z)
        isnothing(result) && return NaN, NaN, NaN

        i_dual,j_dual,k_dual = result

        Ax = _GetFieldValue(domain,X(), datavectorFIT,domain.nodes_u ,domain.edges_v_center,domain.edges_w_center,x, y, z, i, j_dual, k_dual, (i, j, k) -> domain.primal_facets_u[j, k])
        Ay = _GetFieldValue(domain,Y(), datavectorFIT,domain.edges_u_center ,domain.nodes_v,domain.edges_w_center,x, y, z, i_dual, j, k_dual, (i, j, k) -> domain.primal_facets_v[i, k])
        Az = _GetFieldValue(domain,Z(), datavectorFIT,domain.edges_u_center ,domain.edges_v_center,domain.nodes_w,x, y, z, i_dual, j_dual, k, (i, j, k) -> domain.primal_facets_w[i, j])
        return Ax, Ay, Az
end

