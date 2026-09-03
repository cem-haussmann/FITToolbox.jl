# FieldReconstruction.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 21/03/2026

module FieldReconstruction

import ..FITDomain
import ..convert_to_meter 
import ..Normals
import ..X
import ..Y
import ..Z

abstract type TopologicalEntity  end
abstract type EdgeEntity   <: TopologicalEntity end
abstract type FacetEntity  <: TopologicalEntity end
abstract type NodalEntity  <: TopologicalEntity end
abstract type VolumeEntity <: TopologicalEntity end


struct PrimalEdge    <: EdgeEntity  end
struct DualEdge      <: EdgeEntity  end
struct PrimalFacet   <: FacetEntity end
struct DualFacet     <: FacetEntity end
struct PrimalNode    <: NodalEntity end
struct DualNode      <: NodalEntity end

#for fields reconstruction still needs to be inserted
struct PrimalVolume  <: NodalEntity end
struct DualVolume    <: NodalEntity end

include("_InterpolateHelpers.jl")
include("InterpolatePrimalEdgeDualFacet.jl")
include("InterpolateDualEdgePrimalFacet.jl")
include("InterpolateNodalQuantities.jl")

export interpolate
export PrimalEdge, DualEdge, PrimalFacet, DualFacet, PrimalNode, DualNode, PrimalVolume, DualVolume

end # module FieldReconstruction