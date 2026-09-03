# FieldReconstruction.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: long time ago

function create_circular_loop_source(config, radius, center_u, center_v, center_w,
                                     normal::Normals; units="m")
    _create_circular_loop_source(config,radius,center_u,center_v,center_w, normal;units)
end

function _create_circular_loop_source(config,radius,center_u,center_v,center_w, ::Y;units="m")
    throw(ArgumentError("only X() normals are implemented; got Y()"))
end

function _create_circular_loop_source(config,radius,center_u,center_v,center_w, ::Z;units="m")
    throw(ArgumentError("only X() normals are implemented; got Z()"))
end


function _create_circular_loop_source(config, radius, center_u, center_v, center_w, ::X;
                                      units="m")
    Nodes_V, Nodes_W = config.nodes_v, config.nodes_w
    Ev_c, Ew_c = config.edges_v_center, config.edges_w_center
    Nu, Nv, Nw, Np = config.Nu, config.Nv, config.Nw, config.Np
 
    s = check_units(units)
    cu, cv, cw = center_u*s, center_v*s, center_w*s
    R = radius*s
 
    # The loop lies in the v–w plane, so only v and w constrain the radius; u only
    # has to fall inside the domain.
    config.nodes_u[1] <= cu <= config.nodes_u[end-2] ||
        (println("Coil plane lies outside the domain in u-direction"); return NaN)
    cv - R >= Nodes_V[1] && cv + R <= Nodes_V[end-2] ||
        (println("Coil does not fit in domain in v-direction"); return NaN)
    cw - R >= Nodes_W[1] && cw + R <= Nodes_W[end-2] ||
        (println("Coil does not fit in domain in w-direction"); return NaN)
 

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
