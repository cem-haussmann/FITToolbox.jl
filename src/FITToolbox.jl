# FITToolbox.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 21/03/2026

module FITToolbox
using LinearAlgebra
using SparseArrays
using Base.Threads
using TOML

abstract type GridTopology end

"""
    Primal()

Selects the operator on the primal grid, e.g. `get_curl(domain, Primal())` for `C` and
`get_gradient(domain, Primal())` for `G`.
"""
struct Primal <: GridTopology end

"""
    Dual()

Selects the operator on the dual grid, e.g. `get_curl(domain, Dual())` for `C̃ = Cᵀ` and
`get_divergence(domain, Dual())` for `S̃ = -Gᵀ`.
"""
struct Dual   <: GridTopology end

"""
    Direction

An axis of the grid: `DirX()`, `DirY()` or `DirZ()`, along u, v or w. It gives the normal
of a loop source or a cut plane, the axis of a cylinder, torus or helix, the component of
an edge or facet quantity, and, with `Positive()` or `Negative()`, a face of the domain.
"""
abstract type Direction end

"""
    DirX()

The u axis, x in the views; see `Direction`.
"""
struct DirX <: Direction end

"""
    DirY()

The v axis, y in the views; see `Direction`.
"""
struct DirY <: Direction end

"""
    DirZ()

The w axis, z in the views; see `Direction`.
"""
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
include("CreateTorus.jl")
include("CreateHelix.jl")

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
include("CreateCircularLoop.jl")
include("CreateRectangularLoop.jl")
include("CreatePolylineSource.jl")
include("CreateElectricCharge.jl")
include("SourceLengths.jl")

# --- Save domain and load domain ---
include("SaveLoad.jl")

# --- Plotting (implemented in ext/FITToolboxMakieExt.jl) ---
"""
    plot_nodal_values(config, ::Primal, data, normal::Direction; kwargs...)

Plot nodal values on a cut plane through `config` with the given `normal`, as a heatmap
with a colour bar; returns the `Figure`. `data` has one value per primal node (length
`Np`), e.g. a potential `Φ`.

Keywords:
- `pos = 0`, `units = "m"`: the position of the plane along `normal`;
- `label = "Potential Φ [V]"`: the label of the colour bar;
- `clip_range = nothing`: a fixed colour range `(lo, hi)`, else that of the plane;
- `logscale = false`: a log colour scale; zero and negative values are left out, and it
  throws if the plane has no positive value;
- `plot_negative_gradient = false`: field lines of `-G Φ` on top, the electric field of
  a potential, with the `streamplot` settings `density = 1.2`, `stepsize = 0.5`,
  `maxsteps = 1000` and `gridsize = (32, 32)`.

Throws a `DimensionMismatch` if `data` does not have length `Np`.

Requires a Makie backend. Run `using CairoMakie` (or GLMakie / WGLMakie)
before calling; the implementation lives in a package extension.
"""
function plot_nodal_values end

"""
    view_domain(domain, fields...; size = (1400, 750), refine = 1) -> Figure

Interactive 3D view of the objects in `domain`, read from its history: bricks as
boxes, spheres, cylinders, tori and helices as they are, loop sources as their outline and point charges
as dots (red positive, blue negative), in their colours, inside the outline of the
domain. In the voxel view a charge sits on the node that receives it.

The panel next to the scene has
- a view menu with three views:
  - "defined": each object exactly as defined, on its own. A cut-out, a solid created
    with `color = "background"`, is drawn transparent;
  - "combined": the solids put together as smooth surfaces. Where they overlap the
    later one wins, as in the material, and cut-outs are subtracted (and still drawn
    transparent over their holes). Solids of the same colour and material that touch,
    e.g. the pieces of a wire, form one surface without a seam between them. The
    surfaces come from the objects' shapes, not from the FIT grid, sampled at least
    240 times along the longest side of the domain, and finely enough that the
    thinnest solid, e.g. a thin plate, is three samples thick, as far as a limit on
    the number of samples allows;
  - "voxels": the cells each object occupies on the FIT grid, by the same rules: what
    the solver sees.
  Values written directly into `domain.material` are not shown;
- the size of the grid: its number of cells (see `get_cell_count`) and the cells along
  each axis;
- a checkbox per object, in its colour, to show or hide it, or a checkbox per object
  type when there are more than 15 objects, and then a menu of all objects by name,
  with a search field, in which picking an object shows or hides it;
- "Show all"; a click on an object (without dragging) hides it and names it;
- "Reset view", which brings the camera back to where it started;
- a cut: a plane along x, y or z, moved over the grid nodes with a slider, which
  hides everything beyond it. In the voxel view the cut face is filled with the
  cells, showing the staircase in cross-section. The grid lines in the cut plane can
  be shown;
- for loop and polyline sources: the grid edges they occupy, as `get_source` uses
  them, and optionally arrows for the direction of the current;
- with `fields`, a field on the plane of the cut, chosen from a menu.

Each field is given as `"name" => (data, entity)`, where `entity` says what `data` is:
`PrimalNode()` or `DualNode()` for a scalar (length `Np`), `PrimalEdge()`, `DualEdge()`,
`PrimalFacet()` or `DualFacet()` for a vector quantity (length `3Np`), e.g. `ê`, `ĥ`,
`b̂̂` or `d̂̂`. The field is reconstructed with `interpolate` at the grid nodes of the
plane; `refine = 2, 4, …` adds points between them. The plane uses the axis and the
position of the cut, but has its own checkbox, so it can be shown without cutting the
objects. For a vector the plane shows its magnitude or one component; for complex data
the instantaneous value Re(F e^{jφ}) with a phase slider, or its real part, imaginary
part, magnitude or phase. The colours span the whole field by default, so that planes
can be compared, or only the plane ("fit colours to plane"): symmetric about zero for
values of both signs (and for the instantaneous value of complex data), from zero for
magnitudes, and from the smallest to the largest value otherwise, e.g. 0 … 1 V for a
potential between plates at 0 V and 1 V. Its limits can be set: with the two handles of
a slider over that range, in decades on a log scale, or typed into "min" and "max",
also beyond it, e.g. to bring out weak parts of the field or to compare planes at the
same scale. An end left empty, or at the end of the slider, follows the field; "auto"
resets both, as does choosing another field or quantity. The
colour map is `jet`, or `viridis`; values without negative ones, such as magnitudes or
that potential, can be shown on a log scale. The name labels the
colour bar, so it can carry the unit, e.g. `"B (T)"`. Vector fields can also be shown as
arrows: their part in the plane, on a coarser lattice (about 22 along the longer side),
all the same length or scaled to the magnitude, for complex fields at the phase of the
slider.

A field that changes over time, e.g. from FDTD, is given as
`"name" => (snapshots, entity, times)`: `snapshots` is a vector of data vectors, or a
matrix with one per column, recorded at `times` (in s). It can also be any other
`AbstractVector` of data vectors, e.g. one that reads each snapshot from a file when it
is indexed: then only the snapshot shown is held in memory, and the colour range reads
each snapshot once, the first time the field is shown. A slider steps through them and
"play" runs through them, at a speed chosen from 0.5 to 10 snapshots per second
(2 by default); the colour bar label shows the time. Switching to another
time series keeps the time, at its nearest snapshot, so fields recorded half a step
apart (`ê` and `ĥ`) can be compared. The colours span all snapshots by default.

```julia
view_domain(domain, "Φ (V)" => (Φ, PrimalNode()), "B (T)" => (b, PrimalFacet()))
view_domain(domain, "E (V/m)" => (e_snapshots, PrimalEdge(), e_times),
                    "H (A/m)" => (h_snapshots, DualEdge(), h_times))
```

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
export create_domain, get_cell_count
export list_objects
export remove_object!, undo!, change_resolution, get_source

# Operators
export get_curl
export get_gradient
export get_divergence

# Create Objects
export create_sphere!
export create_cylinder!
export create_torus!
export create_helix!, get_helix_ends
#export create_cube! # old call, deprecated, exported by Base.@deprecated
export create_brick!
export create_circular_loop_source!
export create_rectangular_loop_source!
export create_polyline_source!, get_source_ends, get_source_path
export get_source_length
export create_electric_charge!
#export create_circular_loop_source # old call, deprecated, exported by Base.@deprecated

#Option to Load and Save Domain
export save_domain, load_domain

# Material matrices
export get_conductivity
export get_permittivity
export get_reluctivity

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