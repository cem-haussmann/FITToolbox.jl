# FITToolbox.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 21/03/2026

module FITToolbox
#__precompile__(false)
using LinearAlgebra
using SparseArrays
using Base.Threads
using TOML

abstract type GridTopology end

struct Primal <: GridTopology end
struct Dual   <: GridTopology end

abstract type Direction end

struct DirX <: Direction end
struct DirY <: Direction end
struct DirZ <: Direction end


const ε₀ = 8.8541878188e-12   # F/m
const μ₀ = 4π * 1e-7          # H/m

include("CheckUnits.jl")


# --- Domain ---
include("CreateDomain.jl")
include("History.jl")
include("CreateBrick.jl")
include("CreateSphere.jl")
include("CreateCylinder.jl")

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

# --- Save domain and load domain ---
include("SaveLoad.jl")

# --- Plotting (implemented in ext/FITToolboxMakieExt.jl) ---
"""
    plot_nodal_values(config, ::Primal, data, normal::Direction; kwargs...)

Plot nodal values on a cut plane through `config` with the given `normal`.

Requires a Makie backend. Run `using CairoMakie` (or GLMakie / WGLMakie)
before calling; the implementation lives in a package extension.
"""
function plot_nodal_values end

"""
    view_domain(domain; size = (1100, 700)) -> Figure

Interactive 3D view of the objects in `domain`, read from its history: bricks as
boxes, spheres and cylinders as they are, and loop sources as circles, in their
colours, inside the outline of the domain.

The panel next to the scene has
- a view menu: the objects as modelled, or the voxels, i.e. the cells each object
  occupies on the grid (the later object wins where they overlap, as in the
  material). Values written directly into `domain.material` are not shown;
- a checkbox per object, in its colour, to show or hide it, or a checkbox per object
  type when there are more than 15 objects;
- "Show all"; a click on an object (without dragging) hides it and names it;
- "Reset view", which brings the camera back to where it started;
- a cut: a plane along x, y or z, moved over the grid nodes with a slider, which
  hides everything beyond it. In the voxel view the cut face is filled with the
  cells, showing the staircase in cross-section. The grid lines in the cut plane can
  be shown;
- for loop sources: the grid edges they occupy, as `get_source` uses them, and
  optionally arrows for the direction of the current.

Requires a Makie backend with a depth buffer: run `using GLMakie` (window) or
`using WGLMakie` (notebook) before calling. CairoMakie cannot draw overlapping 3D
objects correctly, so `view_domain` throws an error there. The implementation lives
in a package extension.
"""
function view_domain end

# -------------------------------------------------------
# Exports
# -------------------------------------------------------

# Domain
export FITDomain
export create_domain
export list_objects
export remove_object!, undo!, change_resolution, get_source

# Operators
export get_curl
export get_gradient
export get_divergence

# Create Objects
export create_sphere!
export create_cylinder!
#export create_cube! # old call, deprecated, exported by Base.@deprecated
export create_brick!
export create_circular_loop_source!
export create_circular_loop_source   # deprecated

#Option to Load and Save Domain
export save_domain, load_domain

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
export plot_nodal_values, view_domain

export Direction, DirX, DirY, DirZ

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