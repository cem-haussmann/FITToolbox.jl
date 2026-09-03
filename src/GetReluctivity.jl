#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

function get_reluctivity(config)
    #@warn "The create of material matrices was restructured and simplified and still needs to be tested thoroughly"
    
    Edge_U = config.edges_u
    Edge_V = config.edges_v
    Edge_W = config.edges_w

    Nu = config.Nu
    Nv = config.Nv
    Nw = config.Nw

    Material_distribution = config.material

    Mu = 1 
    Mv = Nu
    Mw = Nu*Nv
                
    Np = Nu*Nv*Nw

    #extend a virtual layer to get rid of nasty if statements
    Extended_Edge_U = vcat(Edge_U[1], Edge_U, Edge_U[end]) 
    Extended_Edge_V = vcat(Edge_V[1], Edge_V, Edge_V[end])
    Extended_Edge_W = vcat(Edge_W[1], Edge_W, Edge_W[end])

    #index one corresponds to kappa/sigma
    #index two corresponds to eps_r
    #index three corresponds to mu_r
    function _check_outer_filling(m)
        v = @view m[:, :, :, 3]
        touching = any(!=(1.0), @view(v[1,   :,   :])) || any(!=(1.0), @view(v[end, :,   :])) ||
                   any(!=(1.0), @view(v[:,   1,   :])) || any(!=(1.0), @view(v[:, end,   :])) ||
                   any(!=(1.0), @view(v[:,   :,   1])) || any(!=(1.0), @view(v[:,   :, end]))
        touching && @warn "outermost cells hold non-default permeability; \
                           they are averaged against vacuum at the boundary" maxlog=1
    end
    _check_outer_filling(Material_distribution)

    materialPropertyOfEachVoxel = fill(1.0/μ₀, Nu+1, Nv+1, Nw+1)
    @views materialPropertyOfEachVoxel[2:end-1, 2:end-1, 2:end-1] .*= 1.0 ./ Material_distribution[:, :, :, 3]

    newAveragedMaterialProperty = zeros(Nu*Nv*Nw*3)

    #define the facet areas for weighting of the permittivity
    A_u = (Extended_Edge_V .* Extended_Edge_W' )
    A_v = (Extended_Edge_U .* Extended_Edge_W' )
    A_w = (Extended_Edge_U .* Extended_Edge_V' )

    @threads for k in 1:Nw
        @inbounds for j in 1:Nv, i in 1:Nu
                p = 1 + (i - 1) * Mu + (j - 1) * Mv + (k - 1) * Mw
                
                #u-dir
                newAveragedMaterialProperty[p] = (Extended_Edge_U[i+1]*0.5*materialPropertyOfEachVoxel[i+1,j+1,k+1] + 
                                                    Extended_Edge_U[i]*0.5*materialPropertyOfEachVoxel[i,j+1,k+1]) / A_u[j+1,k+1]        
                #v-dir
                newAveragedMaterialProperty[p+Nu*Nv*Nw] = (Extended_Edge_V[j+1]*0.5*materialPropertyOfEachVoxel[i+1,j+1,k+1] + 
                                                            Extended_Edge_V[j]*0.5*materialPropertyOfEachVoxel[i+1,j,k+1]) / A_v[i+1,k+1]
                #w-dir
                newAveragedMaterialProperty[p+2*Nu*Nv*Nw] = (Extended_Edge_W[k+1]*0.5*materialPropertyOfEachVoxel[i+1,j+1,k+1] + 
                                                                Extended_Edge_W[k]*0.5*materialPropertyOfEachVoxel[i+1,j+1,k]) / A_w[i+1,j+1]
        end
    end
    return spdiagm(0=>newAveragedMaterialProperty)  
end

function M_ν(config)
    return get_reluctivity(config)
end