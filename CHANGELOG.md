# Changelog

All notable changes to FITToolbox.jl are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). As long as
the major version is 0, a change in the minor version marks a breaking release.

## [0.3.0] - Unreleased

This release records every object placed in a domain in a history. Objects can be named,
listed, removed and undone; the whole model can be moved to another grid; and a domain
can be saved to a file and loaded again. The functions that create objects now return
the id of the new object.

### Changed

- **Breaking:** `create_cube!` is renamed to `create_brick!`. The old name still works,
  but is deprecated.
- **Breaking:** `create_brick!` and `create_sphere!` return the `obj_id` of the new
  object. Previously `create_cube!` returned `nothing` and `create_sphere!` the domain.
- **Breaking:** `create_circular_loop_source` is replaced by
  `create_circular_loop_source!`, which records the loop in the domain and returns its
  `obj_id`. The source vector comes from `get_source(domain, obj_id)` or
  `get_source(domain, name)`, and the current is set with the new keyword `current`
  (default 1 A) instead of by scaling the returned vector. The old name still returns
  the vector, but is deprecated.
- Deprecation warnings appear only when Julia runs with `--depwarn=yes`, as in
  `Pkg.test`; by default the old names work silently.

### Added

- **History:** `create_brick!`, `create_sphere!` and `create_circular_loop_source!` are
  recorded in the domain, each with a unique `obj_id`. The material is rebuilt from the
  history whenever it changes.
  - `list_objects(domain)` shows the existing objects as a table of id, type, name,
    geometry and material, and can be indexed like a vector.
  - `remove_object!(domain, obj_id)` and `remove_object!(domain, name)` remove an object.
  - `undo!(domain)` reverts the last step: a created object is removed, a removed one
    comes back.
  - `change_resolution(domain, resolution)` and
    `change_resolution(domain, Edges_U, Edges_V, Edges_W)` return a new domain with the
    same extent and background material and the whole history replayed on the new grid,
    equidistant or not. The original domain is not changed.
  - `get_source(domain, obj_id)` and `get_source(domain, name)` return the source vector
    of one source, discretized on the current grid.
- **Names:** every create function takes a keyword `name`, which must be unique among
  the existing objects. Unnamed objects are called `obj<id>`, a form reserved for them.
- **Colours:** `create_brick!` and `create_sphere!` take a keyword `color`: a hex code
  such as `"#B87333"` or a named colour, including `"copper"`, `"aluminium"` and
  `"gold"`. Colours are stored with the object and shown by `list_objects`; the plotting
  functions do not use them yet.
- **Save and load:** `save_domain(path, domain)` writes the grid, the background
  material and the full history to a human-readable TOML file; `load_domain(path)`
  rebuilds the domain from it, with `undo!` still working. Files are checked on loading:
  ids and names must be unique and every removal must refer to an existing object.
  Values written directly into `domain.material` are not saved, and `save_domain` warns
  if there are any.
- **Complex currents:** `current` may be complex, e.g. a phasor in the frequency domain.
  `get_source` then returns a `ComplexF64` vector; real currents keep `Float64`.
- `create_circular_loop_source!` throws an `ArgumentError` for a radius that is not
  positive, as `create_sphere!` does.
- New dependency: the `TOML` standard library.
- Test files:
  - `test/history.jl`: recording, ids, names, removal, undo, change of resolution,
    `list_objects` and `get_source`.
  - `test/saveload.jl`: round trips, hand-written files and invalid or inconsistent
    files.
  - `test/objects.jl`: material filled by `create_brick!` and `create_sphere!`.
  - `test/sources.jl`: complex currents and the deprecated `create_circular_loop_source`.
  - `test/material_matrices.jl`: the aliases `M_σ`, `M_ε` and `M_ν`.

### Fixed

- The docstrings of `get_index_entity`, `get_reluctivity` and the dual-facet averaging
  were separated from their functions by a blank line and therefore not attached to
  them; they now show up in the help.

### Documentation

- The coil, eddy-current and SPFD dosimetry notebooks use `create_circular_loop_source!`
  with `current` and `get_source`.
- New notebook `examples/snippets/Domain_and_History.ipynb`.

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

[0.3.0]: https://github.com/cem-haussmann/FITToolbox.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/cem-haussmann/FITToolbox.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/cem-haussmann/FITToolbox.jl/releases/tag/v0.1.0
