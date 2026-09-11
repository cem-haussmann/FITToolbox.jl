# Changelog

All notable changes to FITToolbox.jl are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). As long as
the major version is 0, a change in the minor version marks a breaking release.

## [0.2.0] - 2026-09-11

This release corrects the material matrices at the domain boundary. Interior entries are
unchanged, and results that eliminate every boundary degree of freedom (PEC walls,
grounded boundaries) are unaffected. Solutions whose boundary degrees of freedom are live —
homogeneous Neumann, Robin, magnetic walls — change, mostly near the boundary.

### Changed

- **Breaking:** `get_permittivity` and `get_conductivity` give edges lying in a boundary
  plane their geometrically correct half (face) or quarter (domain edge) dual facet,
  consistent with `config.dual_facets_*`. Previously a ghost layer duplicated the
  outermost cells, doubling these facets and averaging the outermost material against
  vacuum.
- **Breaking:** `get_reluctivity` gives facets lying in a boundary plane their half dual
  edge, consistent with `config.dual_edges_*`. Previously the dual edge was doubled.
- **Breaking:** the dead entries in the last grid plane of each direction are now exactly
  zero in all three material matrices: the dead edges in $\mathbf M_\varepsilon$ and
  $\mathbf M_\sigma$, the dead facets in $\mathbf M_\nu$. The matrices are therefore
  positive semi-definite on the full $3N_p$ space. Code that inverts them directly must
  restrict to real entries or use a diagonal pseudo-inverse ($1/0 := 0$).
- **Breaking:** `create_circular_loop_source` throws an `ArgumentError` when the loop does
  not fit in the domain, instead of printing a message and returning `NaN`.
- Magnetic walls ($\mathbf n \times \vec H = 0$) are the natural boundary condition of
  $\mathbf C^\mathsf{T}\mathbf M_\nu\mathbf C$ and are assembled from the plain curl. The
  examples no longer apply a boundary projection from the left, which had moved the three
  negative walls about half a cell into the domain.
- `get_ghost_matrix` is no longer needed to assemble system matrices from the package's
  material matrices, since their dead entries are zero. It remains useful for cleaning raw
  gradient or curl vectors.
- The dual-facet averaging returns a floating-point element type for integer material
  arrays.

### Fixed

- `get_index_entity` no longer adds a component offset for single-component entities
  (`PrimalNode`, `DualVolume`, `DualNode`, `PrimalVolume`) when a normal is passed; it
  previously returned an index outside the node vector for `DirY()` and `DirZ()`.
- Removed an unreachable `nothing` branch from `get_index_entity`. Positions outside the
  domain clamp to the nearest boundary entity and warn when the distance exceeds `atol`,
  as before; this is now documented.
- The field in the plotting extension test was generated with the wrong index order.

### Added

- `get_reluctivity` validates the relative permeability and throws an `ArgumentError`
  for values that are not positive, including `NaN`.
- The dual-facet averaging validates the material index and throws an `ArgumentError`
  for anything other than conductivity or permittivity.
- Docstrings for `get_reluctivity`, `get_index_entity`, `create_circular_loop_source`
  and the dual-facet averaging.
- Test files:
  - `test/material_matrices.jl`: boundary and dead entries against the domain geometry,
    cell-wise reference implementations, anisotropy, element types, the energy of uniform
    fields, parallel-plate capacitors, and the discrete Gauss law.
  - `test/indexing.jl`: position/index round trips for all entity types.
  - `test/interpolation_edges_facets.jl`: exact reproduction of linear fields for edge
    and facet quantities.
  - `test/sources.jl`: closure, orientation and enclosed area of the loop source.

### Documentation

- Magnetostatic, eddy-current and SPFD dosimetry notebooks describe the magnetic wall as
  the natural condition; the electrostatic notebooks explain the dead edges.
- The README example assembles the Poisson operator without `get_ghost_matrix`.

## [0.1.0] - 2026-09-07

Initial registered release.

[0.2.0]: https://github.com/cem-haussmann/FITToolbox.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/cem-haussmann/FITToolbox.jl/releases/tag/v0.1.0
