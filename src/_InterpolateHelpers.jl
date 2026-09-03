#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 21/03/2026

function _find_dual_index(centers::AbstractVector, val::Real)
    n = length(centers)
    n < 2 && throw(ArgumentError("Need at least two dual nodes, got $n"))

    lo, hi = first(centers), last(centers)
    tol = 8 * eps(max(abs(lo), abs(hi)))   # scale-aware, unlike eps()

    if val < lo - tol
        @warn "Value is left of first dual node ($lo), clamping" maxlog=3
        return 1
    elseif val > hi + tol
        @warn "Value is right of last dual node ($hi), clamping" maxlog=3
        return n - 1
    end

    return clamp(searchsortedlast(centers, val), 1, n - 1)
end

function _find_primal_index(nodes::AbstractVector, val::Real)
    if val < first(nodes) || val > last(nodes)
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
    # _clamp_dual_index already handles 1D bounds checking and clamping internally
    i_dual = _find_dual_index(domain.edges_u_center, x)
    j_dual = _find_dual_index(domain.edges_v_center, y)
    k_dual = _find_dual_index(domain.edges_w_center, z)

    return i_dual, j_dual, k_dual
end

function _GetFieldValue(domain::FITDomain,normals::Normals,edgeQuantities::AbstractVector,x_pos_array,y_pos_array,z_pos_array,x,y,z,i,j,k,normalization_)
    offset = 0
    if normals isa Y
        offset+=domain.Np
    elseif normals isa Z
        offset+=2*domain.Np
    end

    Nu = domain.Nu
    Nv = domain.Nv

    Np = domain.Np
    
    Mu = 1
    Mv = Nu
    Mw = Nu*Nv
    
    x_d = (x-x_pos_array[i])/(x_pos_array[i+1]-x_pos_array[i])  
    y_d = (y-y_pos_array[j])/(y_pos_array[j+1]-y_pos_array[j])
    z_d = (z-z_pos_array[k])/(z_pos_array[k+1]-z_pos_array[k]) 

    p0 = 1 + (i - 1) * Mu + (j - 1) * Mv + (k - 1) * Mw + offset;

    # +1 in u-direction
    p1 = p0 + 1;
    # +1 in v-direction
    p2 = p0 + Mv;
    # +1 in v-direction; +1 in u direction
    p3 = p0 + Mv + 1;
    # +1 in w-direction
    p4 = p0 + Mw;
    # +1 in w-direction; +1 in u-direction
    p5 = p0 + Mw + 1;
    # +1 in w-direction; +1 in v-direction
    p6 = p0 + Mw + Mv;
    # +1 in w-direction; +1 in v-direction; +1 in u-direction
    p7 = p0 + Mw + Mv + 1;
    
    c00 = edgeQuantities[p0]/(normalization_(i,j,k))*(1-x_d) + edgeQuantities[p1]/(normalization_(i+1,j,k))*x_d
    c01 = edgeQuantities[p4]/(normalization_(i,j,k+1))*(1-x_d) + edgeQuantities[p5]/(normalization_(i+1,j,k+1))*x_d
    c10 = edgeQuantities[p2]/(normalization_(i,j+1,k))*(1-x_d) + edgeQuantities[p3]/(normalization_(i+1,j+1,k))*x_d
    c11 = edgeQuantities[p6]/(normalization_(i,j+1,k+1))*(1-x_d) + edgeQuantities[p7]/(normalization_(i+1,j+1,k+1))*x_d

    c0 = c00*(1-y_d) + c10*y_d
    c1 = c01*(1-y_d) + c11*y_d

    c = c0*(1-z_d) + c1*z_d    
    return c
end