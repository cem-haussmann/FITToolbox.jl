#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

# deprecated: the grid is given, get_curl(config, Primal()) or get_curl(config, Dual())
@deprecate get_curl(config) get_curl(config, Primal())


"""
    get_curl(config, Primal()) -> SparseMatrixCSC    # C, 3Np × 3Np
    get_curl(config, Dual())   -> SparseMatrixCSC    # C̃ = Cᵀ, 3Np × 3Np

The discrete curl. `C` takes voltages on the primal edges, e.g. `ê` or a vector potential
`â`, to the circulation around each primal facet: Faraday's law reads
`C ê = -d/dt b̂̂`, and `b̂̂ = C â`. `C̃ = Cᵀ` takes voltages on the dual edges, e.g. `ĥ`, to
the circulation around each dual facet: Ampère's law reads `C̃ ĥ = d/dt d̂̂ + ĵ`. Both
are metric-free, with entries ±1, and `C * get_gradient(config, Primal())` and
`get_divergence(config, Primal()) * C` vanish exactly on any grid.

Rows and columns follow the vector layout: the u, v and w components, each of length
`Np` and numbered by the node `p = i + (j-1)Nu + (k-1)NuNv` they start from. Entries on
the last node plane of a direction are dead (see `get_ghost_indices`).

# Example

```julia
C  = get_curl(domain, Primal())
Mν = get_reluctivity(domain)
K  = get_curl(domain, Dual()) * Mν * C        # the curl–curl operator C̃ M_ν C
```
"""
function get_curl(config, ::Primal)
    Nu = config.Nu
    Nv = config.Nv
    Nw = config.Nw
    Np = Nu*Nv*Nw
 
    Pu,Pv,Pw = GetPuPvPw(config)
    C = [spzeros(Np,Np) -1*Pw Pv;
            Pw spzeros(Np,Np) -1*Pu;
            -1*Pv Pu spzeros(Np,Np)] 
    return C
end

#C̃ = transpose(C)
function get_curl(config, ::Dual)
    C = get_curl(config, Primal())
    return sparse(C')
end
