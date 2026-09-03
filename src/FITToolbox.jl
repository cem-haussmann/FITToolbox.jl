# FITToolbox.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 21/03/2026

module FITToolbox
#__precompile__(false)
using LinearAlgebra
using SparseArrays
using Base.Threads

abstract type GridTopology end

struct Primal <: GridTopology end
struct Dual   <: GridTopology end

abstract type Normals end

struct X <: Normals end
struct Y <: Normals end
struct Z <: Normals end


const ε₀ = 8.8541878188e-12   # F/m
const μ₀ = 4π * 1e-7          # H/m

include("CheckUnits.jl")


# --- Domain ---
include("CreateDomain.jl")
include("CreateCube.jl")
include("CreateSphere.jl")

# --- Operators ---
include("_GetPuPvPw.jl")
include("GetCurl.jl")
include("GetGradient.jl")
include("GetDivergence.jl")

# --- Material matrices ---
include("_GetDualFacetAveragedMaterialProperty.jl")
include("GetConductivity.jl")
include("GetPermittivity.jl")
include("GetReluctivity.jl")

# --- Boundary + indexing ---
include("GetIndicesOfBoundaries.jl")

# --- Field reconstruction (submodule) ---
include("FieldReconstruction.jl")
using .FieldReconstruction

include("GetIndexOfEntity.jl")
include("GetPositionOfIndex.jl")

# --- CFL ---
include("GetVacuumCFLTime.jl")

# --- Wire ---
include("CreateWire.jl")

# --- Plotting ---
include("Plotting.jl")

# -------------------------------------------------------
# Exports
# -------------------------------------------------------

# Domain
export FITDomain
export create_domain

# Operators
export get_curl
export C,C̃
export get_gradient
export get_divergence
export G,G̃,S,S̃

# Create Objects
export create_sphere!
export create_cube!
export create_circular_loop_source

# Material matrices
export get_conductivity
export get_permittivity
export get_reluctivity
#export M_ν,M_ε,M_σ

# Field reconstruction
export interpolate
export PrimalEdge, DualEdge, PrimalFacet, DualFacet, PrimalNode, DualNode, PrimalVolume, DualVolume
export Primal, Dual

#find the correct index or position in domain
export get_index_entity, get_position_of_index

# Utils
export get_vacuum_cfl_time

# Plotting
export plot_nodal_values

export Normals, X, Y, Z

###---these are all exports from the GetIndicesofBoundaries
# --- Boundary Types---
export DomainFace, Positive, Negative
export BoundaryComponent, NormalComponent, TangentialComponent, NodalComponent
export DomainBoundary

# --- Matrix Exports ---
export get_boundary_matrix
export get_ghost_matrix

# --- Index Exports) ---
export get_boundary_indices
export get_custom_boundary_indices
export get_ghost_indices
export get_all_tangential_boundary_indices
export get_all_normal_boundary_indices
###

end # module FITToolbox