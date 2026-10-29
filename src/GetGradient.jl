#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

# deprecated: the grid is given, get_gradient(config, Primal()) or get_gradient(config, Dual())
@deprecate get_gradient(config) get_gradient(config, Primal())

"""
    get_gradient(config, Primal()) -> SparseMatrixCSC    # G, 3Np × Np
    get_gradient(config, Dual())   -> SparseMatrixCSC    # G̃, 3Np × Np

The discrete gradient. `G` takes values on the primal nodes, e.g. a potential `Φ`, to
their differences along the primal edges: `ê = -G Φ` are the edge voltages of an
electrostatic field. `G̃` does the same for values on the dual nodes, the cell centres,
along the dual edges. Both are metric-free, with entries ±1, and
`get_curl(config, Primal()) * G` vanishes exactly on any grid. `G̃ = -Sᵀ` with
`S = get_divergence(config, Primal())`.

The rows are the edges in the vector layout of `get_curl`, u, v and w components of
length `Np` each, the edge leaving node `p` in each direction at `p`. The edges that
would leave the last node plane of a direction are dead (see `get_ghost_indices`): they
difference unrelated nodes. The material matrices are zero on them, so `Gᵀ M_ε G` needs
no projection; `get_ghost_matrix` removes them from `G Φ` itself.

# Example

```julia
G = get_gradient(domain, Primal())
L = G' * get_permittivity(domain) * G        # the discrete Poisson operator Gᵀ M_ε G
```
"""
function get_gradient(config, ::Primal)
    Pu,Pv,Pw = GetPuPvPw(config)
    return [Pu;
            Pv;
            Pw]
end

# G̃ = - transpose(S)
function get_gradient(config, ::Dual)
    S = get_divergence(config, Primal())
    return -sparse(S')
end

#function G(config)
#    return get_gradient(config, Primal())
#end

#function G̃(config)
#    return get_gradient(config, Dual())
#end