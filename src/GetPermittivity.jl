#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

"""
    get_permittivity(config) -> Diagonal    # M_ε, 3Np × 3Np

The permittivity matrix `M_ε`: a diagonal matrix with one entry per primal edge, taking
the edge voltages `ê` to the electric fluxes `d̂̂ = M_ε ê` through the dual facets the
edges pierce. Each entry is the dual facet area times the permittivity `ε₀ ε_r`
averaged over the four cells around the edge, each weighted by its share of the facet,
divided by the edge length. `ε_r` is read from `config.material[:, :, :, 2]`.

An edge lying in a face of the domain gets the half dual facet inside it, an edge along
an edge of the domain the quarter, consistent with `config.dual_facets_*`. The dead
edges on the last node plane of each direction get exactly zero, so the matrix is only
positive semi-definite on the full `3Np` space: restrict it to the real edges, or use a
diagonal pseudo-inverse (1/0 := 0), before inverting it. `get_conductivity` builds
`M_σ` the same way.

The matrix is a `Diagonal`, which stores only the `3Np` entries: a third of the memory
of a `SparseMatrixCSC`. Products with the sparse operators, such as `G' * M_ε * G`,
are sparse; `diag(M_ε)` gives the entries as a vector. Indexing it with index vectors,
`M_ε[idx, idx]`, gives a dense matrix: use `Diagonal(diag(M_ε)[idx])` instead, or
`sparse(M_ε)` where a sparse matrix is needed.

# Example

```julia
G = get_gradient(domain, Primal())
L = G' * get_permittivity(domain) * G        # the discrete Poisson operator Gᵀ M_ε G
```
"""
function get_permittivity(config)
        return _get_dual_facet_averaged_material(config,2)
end

# Shorthand as in the equations, not exported: user code names its matrices M_ε
# (`M_ε = get_permittivity(domain)`), which an exported function of that name would break
function M_ε(config)
    return get_permittivity(config)
end