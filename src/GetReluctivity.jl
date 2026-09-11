#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

"""
    get_reluctivity(config) -> SparseMatrixCSC
 
Build the FIT reluctivity matrix: a diagonal matrix with one entry per primal facet,
mapping facet fluxes to the magnetic voltages along the corresponding dual edges.
 
Each entry is the dual edge through that facet, weighted by the reluctivity
ν = 1/(μ₀ μ_r) of the two cells it passes through (half a cell length in each),
divided by the primal facet area. Permeability is read from
`config.material[:, :, :, 3]` and must be positive in every cell; zero, negative
or `NaN` values throw an `ArgumentError`.
 
This is the counterpart of the dual-facet averaging used by
[`get_permittivity`](@ref) and [`get_conductivity`](@ref): there the material is
averaged over a dual facet, here along a dual edge, which is the consistent choice
for a quantity whose normal flux is continuous across material interfaces.
 
# Boundary cells
 
Cells outside the domain are given zero length, so a facet lying in a boundary
plane receives the geometrically correct half dual edge, consistent with
`config.dual_edges_*`, and no material from outside the domain enters.
 
The dead facets in the last plane of each direction, i.e. the tangential facet
components on the three positive faces, have no dual edge and receive exactly
zero. The matrix is therefore only positive semidefinite on the full 3Nₚ space:
restrict to real facets, or use a diagonal pseudo-inverse (1/0 := 0), before
inverting it.
 
# Example
 
```julia
Mν = get_reluctivity(domain)
A  = get_curl(domain, Dual()) * Mν * get_curl(domain, Primal())   # curl–curl operator
```
"""

function get_reluctivity(config)    
    Edge_U = config.edges_u
    Edge_V = config.edges_v
    Edge_W = config.edges_w

    Nu = config.Nu
    Nv = config.Nv
    Nw = config.Nw

    Mu = 1 
    Mv = Nu
    Mw = Nu*Nv
                
    all(>(0), @view config.material[:, :, :, 3]) ||
    throw(ArgumentError("relative permeability must be positive in every cell"))

    function _extend(m)
        e = zeros(Nu+1, Nv+1, Nw+1)
        @views e[2:end-1, 2:end-1, 2:end-1] .= 1.0 ./ (μ₀ .* m[:, :, :, 3])
        return e
    end

    materialPropertyOfEachVoxel = _extend(config.material)
    newAveragedMaterialProperty = zeros(Nu*Nv*Nw*3)

    # dual edge lengths (numerator): zero-length ghost cells → half dual edge at the boundary
    Extended_Edge_U = vcat(0.0, Edge_U, 0.0) 
    Extended_Edge_V = vcat(0.0, Edge_V, 0.0)
    Extended_Edge_W = vcat(0.0, Edge_W, 0.0)

    # facet areas (divisors): duplicated end widths keep dead facets at 0/A instead of 0/0
    Extended_Edge_U_A = vcat(Edge_U[1], Edge_U, Edge_U[end])
    Extended_Edge_V_A = vcat(Edge_V[1], Edge_V, Edge_V[end])
    Extended_Edge_W_A = vcat(Edge_W[1], Edge_W, Edge_W[end])

    A_u = Extended_Edge_V_A .* Extended_Edge_W_A'
    A_v = Extended_Edge_U_A .* Extended_Edge_W_A'
    A_w = Extended_Edge_U_A .* Extended_Edge_V_A'

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