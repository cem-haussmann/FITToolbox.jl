#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 06/08/2026

"""
    get_divergence(config, Primal()) -> SparseMatrixCSC    # S, Np × 3Np
    get_divergence(config, Dual())   -> SparseMatrixCSC    # S̃, Np × 3Np

The discrete divergence. `S` takes the fluxes through the primal facets, e.g. `b̂̂`, to the
net outflow of each primal cell; `S̃` takes the fluxes through the dual facets, e.g. `d̂̂`
or a current `ĵ`, to the net outflow of each dual cell, i.e. of each primal node. Both are
metric-free, with entries ±1, and `S * get_curl(config, Primal())` and
`S̃ * get_curl(config, Dual())` vanish exactly on any grid. `S̃ = -Gᵀ` with
`G = get_gradient(config, Primal())`: the discrete Gauss law is `S̃ d̂̂ = q`, and a current
loop is closed when `S̃ ĵ = 0`.

The rows are the cells, numbered by the node they start from (`S`) or the node they
surround (`S̃`), `p = i + (j-1)Nu + (k-1)NuNv`; the columns are the facets in the vector
layout of `get_curl`, u, v and w components of length `Np` each. The primal cells
starting on the last node plane of a direction do not exist; their rows are dead.
"""
function get_divergence(config, ::Primal)
    Pu,Pv,Pw = GetPuPvPw(config)
    return [Pu Pv Pw]
end

# S̃ = - transpose(G)
function get_divergence(config, ::Dual)
    G = get_gradient(config, Primal())
    return -sparse(G')
end

# S̃J of a current kept sparse, without S̃ or a vector of the whole grid: J as pairs
# edge => value, S̃J as node => value at the nodes those edges touch. The net outflow of
# the dual cell around each node is the edge leaving it minus the edge arriving at it in
# each direction, (S̃J)_p = Σ_d J_d[p] - J_d[p - M_d]: an edge takes its current from its
# lower node to its upper one, M_d further on. The sources check with it that they are
# closed; its work and memory grow with the number of edges, not with the grid.
function _dual_divergence(config, J)
    Np = config.Np
    M = (1, config.Nu, config.Nu * config.Nv)
    div = Dict{Int,Float64}()
    for (edge, value) in J
        lower = mod1(edge, Np)
        upper = lower + M[(edge - 1) ÷ Np + 1]
        div[lower] = get(div, lower, 0.0) + value
        upper <= Np && (div[upper] = get(div, upper, 0.0) - value)
    end
    return div
end

#function S(config)
#    return get_divergence(config, Primal())
#end

#function S̃(config)
#    return get_divergence(config, Dual())
#end