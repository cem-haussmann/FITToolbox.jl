#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

"""
    get_conductivity(config) -> Diagonal    # M_σ, 3Np × 3Np

The conductivity matrix `M_σ`: a diagonal matrix with one entry per primal edge, taking
the edge voltages `ê` to the conduction currents `ĵ = M_σ ê` through the dual facets.
It is built like `get_permittivity`, from `σ` in `config.material[:, :, :, 1]` and
without `ε₀`: the dual facet area times the conductivity of the four cells around the
edge, each weighted by its share of the facet, divided by the edge length. Edges with
only non-conducting cells around them get zero, as do the dead edges on the last node
plane of each direction. Like `get_permittivity`, it is a `Diagonal`.

# Example

```julia
# time-harmonic eddy currents: (C̃ M_ν C + jω M_σ) â = ĵ
A = get_curl(domain, Dual()) * get_reluctivity(domain) * get_curl(domain, Primal()) +
    im * 2π * 50 * get_conductivity(domain)
```
"""
function get_conductivity(config)
    return _get_dual_facet_averaged_material(config,1)
end

# Shorthand as in the equations, not exported: user code names its matrices M_σ
# (`M_σ = get_conductivity(domain)`), which an exported function of that name would break
function M_σ(config)
    return get_conductivity(config)
end