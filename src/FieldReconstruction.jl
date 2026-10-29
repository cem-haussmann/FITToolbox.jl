# FieldReconstruction.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 21/03/2026

module FieldReconstruction

import ..FITDomain
import ..convert_to_meter 
import ..Direction
import ..DirX
import ..DirY
import ..DirZ
import .._within

abstract type TopologicalEntity  end
abstract type EdgeEntity   <: TopologicalEntity end
abstract type FacetEntity  <: TopologicalEntity end
abstract type NodalEntity  <: TopologicalEntity end
abstract type VolumeEntity <: TopologicalEntity end


"""
    PrimalEdge()

The edges of the primal grid, from node to node, in vectors of length `3Np` (see
`get_curl` for the layout): the home of the electric voltages `ê` (V) and of the
integrated vector potential `â` (Vs). `interpolate` divides by the edge length, giving
e.g. `E` in V/m.
"""
struct PrimalEdge    <: EdgeEntity  end

"""
    DualEdge()

The edges of the dual grid, from cell centre to cell centre through a primal facet, in
vectors of length `3Np`: the home of the magnetic voltages `ĥ` (A). `interpolate` divides
by the dual edge length, giving `H` in A/m.
"""
struct DualEdge      <: EdgeEntity  end

"""
    PrimalFacet()

The facets of the primal grid, the faces of the cells, in vectors of length `3Np`: the
home of the magnetic fluxes `b̂̂` (Vs). `interpolate` divides by the facet area, giving `B`
in T.
"""
struct PrimalFacet   <: FacetEntity end

"""
    DualFacet()

The facets of the dual grid, each pierced by a primal edge, in vectors of length `3Np`:
the home of the electric fluxes `d̂̂` (C) and of the currents `ĵ` (A). `interpolate`
divides by the facet area, giving `D` in C/m² or `J` in A/m².
"""
struct DualFacet     <: FacetEntity end

"""
    PrimalNode()

The nodes of the primal grid, in vectors of length `Np`, numbered
`p = i + (j-1)Nu + (k-1)NuNv`: the home of the electric potential `Φ` (V).
"""
struct PrimalNode    <: NodalEntity end

"""
    DualNode()

The nodes of the dual grid, the centres of the primal cells, in vectors of length `Np`
numbered by the cell's first node; the entries on the last node plane of each direction
have no cell and are dead.
"""
struct DualNode      <: NodalEntity end

"""
    PrimalVolume()

The cells of the primal grid, in vectors of length `Np` numbered like `DualNode()`: the
home of quantities integrated over the cells. `interpolate` divides by the cell volume,
giving a density.
"""
struct PrimalVolume  <: VolumeEntity end

"""
    DualVolume()

The cells of the dual grid, one around each primal node, in vectors of length `Np`: the
home of the charges `q` (C) of the Poisson problem `Gᵀ M_ε G Φ = q`. `interpolate`
divides by the dual cell volume, giving the charge density in C/m³.
"""
struct DualVolume    <: VolumeEntity end

"""
    interpolate(domain, entity, data, x, y, z; units = "m")

The field described by the FIT vector `data` at the point `(x, y, z)`, given in `units`.
`entity` says where `data` lives:

- `PrimalNode()` or `DualNode()`: values on the nodes, e.g. a potential `Φ`; returns a
  number.
- `PrimalEdge()`, `DualEdge()`, `PrimalFacet()` or `DualFacet()`: values integrated over
  the entities, e.g. `ê`, `ĥ`, `b̂̂` or `d̂̂`; returns the field vector `(Fu, Fv, Fw)`. Each
  component is divided by the length or area of its entity and interpolated between
  its own staggered positions.
- `PrimalVolume()` or `DualVolume()`: values integrated over cells, e.g. a charge `q`;
  returns the density, the value divided by the cell volume.

`data` holds one value per entity, `Np` for nodes and volumes and `3Np` for edges and
facets, and may be complex. The interpolation is trilinear, so a linear field is
reproduced exactly. Within half a cell of the boundary a component may have its
positions on one side only; it is then extrapolated from the outermost pair. Along an
axis with a single cell the quantities centred in the cell are constant.

Returns `NaN`, or a tuple of `NaN`, with a warning for a point outside the domain.
Throws a `DimensionMismatch` if the length of `data` does not fit `entity`.

# Example

```julia
interpolate(domain, PrimalNode(), Φ, 0.5, 0.5, 0.25)                    # the potential, V
Eu, Ev, Ew = interpolate(domain, PrimalEdge(), -G * Φ, 50, 50, 25; units = "cm")   # V/m
```
"""
function interpolate end

include("_InterpolateHelpers.jl")
include("InterpolatePrimalEdgeDualFacet.jl")
include("InterpolateDualEdgePrimalFacet.jl")
include("InterpolateNodalQuantities.jl")

export interpolate
export PrimalEdge, DualEdge, PrimalFacet, DualFacet, PrimalNode, DualNode, PrimalVolume, DualVolume

end # module FieldReconstruction