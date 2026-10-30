# Changelog

All notable changes to FITToolbox.jl are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). As long as
the major version is 0, a change in the minor version marks a breaking release.

## [0.3.0] - Unreleased

This release records every object placed in a domain in a history. Objects can be named,
listed, removed and undone; the whole model can be moved to another grid; and a domain
can be saved to a file and loaded again. The functions that create objects now return
the id of the new object. New are cylinders, tori, helices, loop sources along every axis, and an
interactive 3D view of a model, as modelled and as the grid sees it.

### Changed

- **Breaking:** `create_cube!` is renamed to `create_brick!`. The old name still works,
  but is deprecated.
- **Breaking:** `create_brick!` and `create_sphere!` return the `obj_id` of the new
  object. Previously `create_cube!` returned `nothing` and `create_sphere!` the domain.
- **Breaking:** `create_circular_loop_source` is replaced by
  `create_circular_loop_source!`, which records the loop in the domain and returns its
  `obj_id`. The source vector comes from `get_source(domain, obj_id)` or
  `get_source(domain, name)`, and the current is set with the new keyword `current`
  (default 1 A) instead of by scaling the returned vector. Its arguments follow the
  other create functions, the centre first and then the radius:
  `create_circular_loop_source!(domain, u, v, w, radius, normal)`. The old name still
  returns the vector and takes the radius first, as before, but is deprecated.
- **Breaking:** `get_permittivity`, `get_conductivity` and `get_reluctivity` return a
  `Diagonal` instead of a `SparseMatrixCSC`. It stores only the `3Np` entries, a third
  of the memory, and products with the operators, such as `G' * M_ε * G` and
  `C' * M_ν * C + im * ω * M_σ`, are the same sparse matrices as before and as fast or
  faster. `diag(M)`, `M * x` and `M[p, p]` work as before. Code that indexes a material
  matrix with index vectors, `M[idx, idx]`, now gets a dense matrix: use
  `Diagonal(diag(M)[idx])`, or `sparse(M)` where a `SparseMatrixCSC` is needed.
- Deprecation warnings appear only when Julia runs with `--depwarn=yes`, as in
  `Pkg.test`; by default the old names work silently. This now holds for every old
  name: the one-argument `get_curl(config)` and `get_gradient(config)` are deprecated in
  favour of `get_curl(config, Primal())` and `get_gradient(config, Primal())` instead
  of warning on every call.
- `create_brick!` clips a brick that reaches beyond the domain on any side, with a
  warning, as `create_sphere!` and the other solids do, instead of throwing; its corner
  may lie outside, e.g. for a plate through the whole domain. It still throws if a length
  is not positive and finite.
- A brick that fills no cell, because less than about half a cell of it lies in the
  domain along an axis (e.g. a foil thinner than half a cell), is recorded with the
  warning "no cells were filled", as a sphere, cylinder, torus or helix thinner than a
  cell is, instead of throwing. It fills its cells on a finer grid, so `change_resolution`
  to a coarser grid and back no longer fails on it.
- The centre of a sphere or torus and the base of a cylinder or helix may lie outside the
  domain, as the corner of a brick may, e.g. for a sphere centred just beyond a wall: the
  solid is clipped to the part that reaches in, with a warning. It used to be an error
  even when most of the solid lay inside, while a sphere centred on the wall was fine.
  Every solid that misses the domain altogether throws "the sphere lies outside the
  domain: it spans u = 4.7 … 5.3 m, the domain u = 0.0 … 1.0 m", also a brick, which
  used to report a zero extent after snapping instead.
- A solid that does not fit the grid when `change_resolution` or `load_domain` replays the
  history names itself, as the sources do: "object "ball" (obj_id 2) does not fit this
  grid: the sphere lies outside the domain: …", and `load_domain` puts the file in
  front. It used to say only "the sphere lies outside the domain", whichever sphere it was.
- **Breaking:** materials are checked: the solids and the background of `create_domain`
  throw an `ArgumentError` for a conductivity `σ` that is negative and for `ε_r` or `μ_r`
  that are not positive, and for any of them that is not finite. A zero `μ_r` used to
  surface only in `get_reluctivity`, a negative `σ` or `ε_r` never, but as a solver that
  failed or a lossy medium gaining energy. `load_domain` checks a file the same way.
- Lengths of solids and loops must be finite as well as positive, and currents and
  charges finite: an infinite radius used to fill or clip against the whole domain
  without an error. `load_domain` rejects `inf` and `nan` in a file, naming the place.
- The errors of a loop that does not fit name the range it must lie within and the
  value that was given, e.g. "the coil plane does not fit in the domain in u-direction:
  it must lie within u = 0.0 … 0.8 m (nodes[1] … nodes[end-2]), got u = 0.9 m", instead
  of calling a plane in the last two node planes outside the domain.
- **Breaking:** a face of a brick exactly halfway between two nodes snaps to the lower
  node, as the corners of the loop sources and the points of polylines and charges do.
  It used to snap to the upper one, so a rectangular loop around a brick given by the
  same coordinates could sit a cell off.
- `get_gradient(config, Dual())` and `get_divergence(config, Dual())` return a
  `SparseMatrixCSC`, as `get_curl(config, Dual())` does, instead of a lazy adjoint.
- The keywords `name` and `color` take any string, e.g. a `SubString`.
- The closure checks of the sources are errors instead of warnings: every loop is
  closed and every polyline connected by construction, so a failed check is a bug in
  the package. They compute `S̃J` from the edges of the source alone, without `S̃` or a
  vector of the whole grid.
- Wire sources (circular and rectangular loops, polylines) are kept as the list of edges
  they occupy, whose memory grows with the length of the wire, not with the grid:
  creating one, `get_source_path`, `get_source_length`, the checks of
  `change_resolution` and `load_domain`, and `view_domain` drawing its grid edges build
  no vector of the whole grid, and `get_source` builds only its result. A circular loop
  visits only the cells around the circle, no longer its whole plane. On a grid of
  2.9 M nodes a wire of 15 edges took about 197 MiB in each of these calls and a circle
  of 116 edges 88–153 MiB, now 10–20 KiB each, and `get_source` its 66 MiB result; 100
  circles, e.g. the turns of a coil, took 15 GiB to create, now about 1 MiB. Creating a
  point charge builds no such vector either.
- Spheres, cylinders, tori and helices fill the material by scanning only the cells
  around them instead of the whole grid, with the same result cell for cell: creating
  them, and replaying them in `change_resolution`, `load_domain` and the voxel view of
  `view_domain`, now costs in proportion to their size on the grid. On 8 M nodes, a
  model of five solids took 4.3 s to create, now 0.9 s; opening the voxel view 5.1 s,
  now 2.4 s.
- `create_domain(domain_size, resolution)` throws an `ArgumentError` when either does
  not have three values, instead of warning and returning `NaN`.
- `interpolate` no longer warns for points within half a cell of the boundary. There, a
  dual component has a dual node on one side only and is extrapolated from the
  outermost pair, which is exact for linear fields; the warning called this "clamping".

### Added

- **History:** `create_brick!`, `create_sphere!`, `create_cylinder!`, `create_torus!`,
  `create_helix!`,
  `create_circular_loop_source!`, `create_rectangular_loop_source!`,
  `create_polyline_source!` and `create_electric_charge!` are recorded in the domain,
  each with a unique `obj_id`. The material is rebuilt from the
  history whenever it changes.
  - `list_objects(domain)` shows the existing objects as a table of id, type, name,
    geometry and material, and can be indexed like a vector.
  - `remove_object!(domain, obj_id)` and `remove_object!(domain, name)` remove an object.
  - `undo!(domain)` reverts the last step: a created object is removed, a removed one
    comes back.
  - `change_resolution(domain, resolution)` and
    `change_resolution(domain, Edges_U, Edges_V, Edges_W)` return a new domain with the
    same extent and background material and the whole history replayed on the new grid,
    equidistant or not. The original domain is not changed. A resolution must divide
    the extent; otherwise an `ArgumentError` names the direction.
  - `get_source(domain, obj_id)` and `get_source(domain, name)` return the source vector
    of one source, discretized on the current grid.
- **Names:** every create function takes a keyword `name`, which must be unique among
  the existing objects. Unnamed objects are called `obj<id>`, a form reserved for them.
- **Colours:** `create_brick!`, `create_sphere!`, `create_cylinder!`, `create_torus!` and
  `create_helix!` take a keyword
  `color`: a hex code such as `"#B87333"` or a named colour, including `"copper"`,
  `"aluminium"` and `"gold"`. Colours are stored with the object, shown by
  `list_objects` and used by `view_domain`. `color = "background"` makes a solid a
  cut-out: it takes the background material of the domain, e.g. to cut a hole into an
  earlier object, and `view_domain` shows it as a hole.
- **Cylinders:** `create_cylinder!(domain, u_b, v_b, w_b, radius, height, axis)`, a solid
  circular cylinder from the centre of its base face along `DirX()`, `DirY()` or
  `DirZ()`. Cells are filled by the same 50 % rule as for spheres.
- **Tori:** `create_torus!(domain, u_c, v_c, w_c, R, r, axis)`, a solid ring around its
  centre with ring radius `R` and tube radius `r < R`, its symmetry axis along `DirX()`,
  `DirY()` or `DirZ()`. Same 50 % rule; saved and loaded, and drawn by `view_domain`.
- **Helices:** `create_helix!(domain, u_b, v_b, w_b, R, r, pitch, turns, axis)`, a solid
  helical wire, e.g. a winding: radius `R` around the axis, wire radius `r`, rising by
  `pitch` per turn, with round ends so that leads connect. `start_angle` and
  `left_handed` set where it starts and how it turns; `get_helix_ends` returns its end
  points, to attach leads. Warns when the gap between the turns is below a cell, since
  the turns would merge on the grid.
- **Loop sources along every axis:** `create_circular_loop_source!` takes `DirY()` and
  `DirZ()` normals as well as `DirX()`; the current runs counter-clockwise about the
  normal, so the magnetic moment points along it.
- **Rectangular loops:** `create_rectangular_loop_source!(domain, u_o, v_o, w_o, length_1,
  length_2, normal)` (the lengths also as a tuple), a corner and two sides as for
  `create_brick!`. The sides run along grid lines, so there is no staircase: the loop is
  exactly the snapped rectangle, closed, with the moment exactly current × area along the
  normal. Both loop types work alike in `get_source`, `get_source_length` and
  `view_domain`.
- **Polyline sources:** `create_polyline_source!(domain, points; closed = false)`, a wire
  along straight pieces between the given points. Each piece becomes a path of primal
  edges through the nodes nearest to the line, connected by construction (one edge per
  step, within half a cell diagonal of the line), so a straight wire is two points. The
  wire is divergence free except at the ends of an open polyline, where the current
  enters and leaves the grid; an info message names them, and
  `get_source_ends(domain, name)` returns them as positions and node indices. A polyline
  whose edges cancel (there and back along one line) is rejected.
  `get_source_path(domain, name)` returns the edges in their order along the wire, with
  the sign of each (along or against its axis), its length, the nodes passed and the
  piece of the polyline it belongs to, e.g. to build a transmission-line model of the
  wire on the grid. Steps back along the edge just walked, which a sharp corner produces
  and which cancel in the source vector, are left out, so the path holds exactly the
  edges of `get_source`; a polyline that runs along an edge twice has no single order and
  is rejected. Its `planned_lengths` give the length of the planned wire each edge stands
  for: a staircase is longer than the line it follows, √2 at 45° in a plane, so per-length
  quantities of the real wire put on the edges with these lengths keep its electrical
  length, piece by piece, also for a wire that is partly straight and partly diagonal.
  Each piece's length is shared among its own edges. Where parts of the wire lie within a
  cell of each other, their edges cancel and their length goes with them, instead of
  stretching the rest of the wire; a warning says so when more than 5 % is lost.
- **Source lengths:** `get_source_length(domain, name)` returns the length of a polyline
  source, rectangular loop (its sides as well) or circular loop, as planned and as on
  the grid (the edges the source occupies). A staircase is the sum
  of the line's extents along the axes, e.g. `8R` for a circle of radius `R`, on any
  grid.
- **Point charges:** `create_electric_charge!(domain, u, v, w; charge)` records a point
  charge in the history; `charge` is required. `get_source` returns its charge
  vector of length `Np`, the right-hand side of the Poisson problem, with the charge in
  the dual volume of the nearest node.
- **Domain viewer:** `view_domain(domain)` shows the model interactively in 3D. It needs
  a Makie backend with a depth buffer (GLMakie or WGLMakie) and refuses CairoMakie.
  - A view menu switches between the objects as defined, the objects combined into
    smooth surfaces (later objects win, cut-outs subtracted), and the voxels, i.e. the
    cells each object occupies on the grid, drawn as surface meshes on the real,
    possibly graded, cell boundaries.
  - In the combined view, solids of the same colour and material that touch, e.g. the
    leads, elbows and risers of a wire, are one surface, without the creases and
    overlapping walls of separate surfaces where they meet. The sampling resolves the
    thinnest solid with three samples, e.g. 1 mm for 3 mm capacitor plates, within a
    limit of 4 M samples around the solids, and 240 along the longest side at least.
  - A checkbox per object, in its colour (per type above 15 objects), hides or shows it;
    so does a click on the object. "Show all" and "Reset view" restore the defaults.
    With a checkbox per type, a menu of all objects by name, type and obj_id, with a
    search field, shows or hides each one and marks the hidden ones.
  - A cut along x, y or z hides everything beyond a grid plane; in the voxel view the
    cut face is filled with the cells. The grid lines in the cut plane can be shown.
  - Loop and polyline sources can be drawn on the grid edges they occupy, with arrows
    for the direction of the current.
  - Cut-outs (`color = "background"`) are holes in the combined and voxel views, and
    drawn transparent where they are.
  - Fields: `view_domain(domain, "Φ (V)" => (Φ, PrimalNode()), "B (T)" => (b, PrimalFacet()))`
    colours a plane with a nodal, edge or facet quantity, reconstructed with
    `interpolate`. The plane has its own checkbox and follows the cut's axis and
    position; a menu picks the field, and for vectors the magnitude or a component.
    Complex fields show the instantaneous value Re(F e^{jφ}) with a phase slider, or
    their real or imaginary part, magnitude or phase. The colours span the whole field
    or only the plane: symmetric about zero for values of both signs, from zero for
    magnitudes, from the smallest to the largest value otherwise. Values without
    negative ones can be shown on a log scale; jet by default, viridis as an option, and
    a slider narrows the range. `refine` samples the plane between the grid nodes.
  - Time series, e.g. from FDTD: `"E (V/m)" => (snapshots, PrimalEdge(), times)`, with a
    time slider and "play"; switching fields keeps the time at the nearest snapshot.
  - Arrows for vector fields: the part in the plane, on a coarser lattice, all the same
    length or scaled to the magnitude; they follow the phase and the time.
  - Options that do not apply are greyed out and ignore clicks, keeping their state:
    the grid without a cut, arrows for scalars, their length without arrows, the log
    scale for values of both signs.
- **Save and load:** `save_domain(path, domain)` writes the grid, the background
  material and the full history to a human-readable TOML file; `load_domain(path)`
  rebuilds the domain from it, with `undo!` still working. Files are checked on loading:
  ids and names must be unique and not empty, every removal must refer to an existing
  object, and a value that is missing or of the wrong kind (a position without three
  coordinates, an unknown axis or colour, a cell width that is not positive) is an
  `ArgumentError` that names the file and the place in it. Values written directly into
  `domain.material` are not saved, and `save_domain` warns if there are any.
- **Complex currents:** `current` may be complex, e.g. a phasor in the frequency domain.
  `get_source` then returns a `ComplexF64` vector; real currents keep `Float64`.
- `interpolate(domain, PrimalVolume(), …)` and `interpolate(domain, DualVolume(), …)`:
  the density of a quantity integrated over the cells, e.g. a charge `q` from
  `get_source` in C/m³. Both entities were exported, but `interpolate` had no method
  for them.
- `create_circular_loop_source!` throws an `ArgumentError` for a radius that is not
  positive, as `create_sphere!` does.
- `get_cell_count(domain)`: the number of cells of the grid. `view_domain` shows it in its
  panel, with the cells along each axis.
- New dependency: the `TOML` standard library.
- Test files:
  - `test/history.jl`: recording, ids, names, removal, undo, change of resolution,
    `list_objects` and `get_source`.
  - `test/saveload.jl`: round trips, hand-written files, and invalid or inconsistent
    files and files with values missing or of the wrong kind.
  - `test/objects.jl`: material filled by `create_brick!`, `create_sphere!`,
    `create_cylinder!`, `create_torus!` and `create_helix!`; bricks clipped at the
    boundary.
  - `test/sources.jl`: circular and rectangular loops (closure, orientation, the three
    normals, the limits of the domain on either side), polylines, point charges, their
    ends, paths and lengths, complex currents, and the deprecated
    `create_circular_loop_source`.
  - `test/indexing.jl`, `test/interpolation_edges_facets.jl`: the outermost half cells,
    volume quantities, domains one cell thick, and the messages for indices without an
    entity.
  - `test/runtests.jl`: the deprecated names, a docstring for every exported name, the
    boundary indices and projections for all components, and `plot_nodal_values` along
    every axis.
  - The CI runs the tests with two threads, so that the threaded loops run in parallel.
  - `test/viewer.jl`: `view_domain` built with CairoMakie and every control operated:
    colours, groups, cut, voxel labels and meshes, loop edges, and the field plane:
    checks of the fields, sampling, complex modes, colour ranges, time series and the
    panel.
  - `test/material_matrices.jl`: the aliases `M_σ`, `M_ε` and `M_ν`.

- A `FITDomain` displays as a summary of its grid, background and number of objects,
  e.g. `FITDomain: 10×8×6 cells (480), 693 nodes, equidistant`, and as
  `FITDomain(10×8×6 cells, 1 object)` inside a vector. The default display printed every
  field, thousands of characters on one line in the REPL.
- `get_boundary_matrix` and `get_custom_boundary_indices` take the faces as a tuple as
  well as a vector, and `create_polyline_source!` its points.
- Tests:
  - the orientation of `C` and `C̃`: Stokes' theorem for a uniform `B` from its vector
    potential and a uniform current from its `H`, exactly, on a graded grid. The
    identities `C G = 0` and `S C = 0` hold for `-C` too;
  - `get_vacuum_cfl_time` against the leapfrog stability limit, from the largest
    eigenvalue of the system: equal to it on an equidistant grid, below it on a graded
    one;
  - the material, length, current and charge checks, the brick clipped below the
    domain, the loop errors, the display of a domain, tuples of faces and points;
  - `interpolate` on complex data, for nodes and for every edge and facet entity;
  - `create_domain` for a size smaller than half a cell and for a resolution that does
    not divide the size, and the size checks of the `FITDomain` constructor. Every error
    and warning of the package is now reached by a test, except the three checks that
    report a bug in the package itself;
  - Aqua.jl: method ambiguities, type piracy, unbound type parameters, stale
    dependencies and compat entries. `[compat]` now covers the test dependencies.
  - `test/manual_glmakie.jl`, run by hand on a machine with a display, not by CI: the
    views of `view_domain` drawn by GLMakie in hidden windows, a click picking the
    object under the mouse, "play" and its speed, the colour limits, and the time per frame with
    1000 objects; a PNG of each view for a look by eye.

- `view_domain` on large models (measured on 8 M cells, and with 1000 objects):
  - With more than 15 objects, objects of one colour are drawn as one plot per view
    instead of one plot each: 13 plots instead of 3004 for 1000 objects, and the
    figure is built in 1.3 s instead of 16 s. Drawing in GLMakie is 25 times faster: a
    frame of 1000 spheres took 220 ms as 3000 plots and takes 9 ms merged. Each object
    is still shown, hidden and named by a click on its own; a hidden one's points are
    NaN in the merged plot, which is not drawn.
  - The voxel view finds the faces of the objects once; moving the cut keeps those
    before it and closes the cut, instead of scanning the whole grid on every step of
    the slider (0.03 s instead of 0.08 s per step, growing with the surface of the
    objects, not with the number of cells).
  - The voxel labels are written directly into an array of the smallest integer type,
    one byte per cell up to 255 objects, instead of replaying the model into a copy of
    the domain: 8 MiB instead of 216 MiB.
  - The colour range of a field covers every grid node of every snapshot. It took every
    s-th node above a million in all, and missed a peak narrower than that, e.g. at a
    source edge, which then saturated the colours.
  - A time series given as a vector that reads its snapshots on demand, e.g. from a
    file, is read on demand: opening the view reads the first snapshot, showing a step
    reads that one. Every snapshot used to be read at once (40 snapshots of 24 MiB:
    984 MiB to open the view, now 64 MiB).
- The limits of the colour range in `view_domain` can be set: with the two handles of a
  slider over the automatic range, in decades on a log scale, or typed into "min" and
  "max", also beyond that range. An end left empty, or at the end of the slider,
  follows the field; "auto" resets both, as does another field, quantity or complex
  mode. The slider that narrowed the range by a factor towards zero is gone; scaled
  arrows follow the upper limit as they followed it.
- "play" in `view_domain` runs at a speed chosen under it, 0.5 to 10 snapshots per
  second, 2 by default, also changed while it plays. It ran at up to 10, too fast for a
  series of a few snapshots.
  Pressing pause and then play within one step no longer stops the new run.

### Fixed

- The voxel view of `view_domain` kept its cell labels as `Int` although they are made as
  `Int32` to halve their memory; they stay `Int32` now.
- A test of the combined view printed the warning of a sphere smaller than a cell; it
  checks the warning now.
- A circular loop on a grid too coarse to resolve it on one side, where its crossing
  edges fall on the line through its centre, was not closed: the direction of an edge
  was taken from its position relative to the centre, which is ambiguous there. It now
  follows from the side of the edge the inside is on, which is the same for resolved
  loops and closes the loop on any grid.
- A circular loop reaching into the first cells of the domain, within half a cell of
  `nodes[1]`, was not closed either: the edges on the first node line were never
  marked, and only a warning said so. Every loop the checks accept is now closed.
- `change_resolution` and `load_domain` check every source on the grid they produce. A
  loop that no longer fits within `nodes[1] … nodes[end-2]`, or a polyline whose points
  all snap to one node, used to pass unnoticed and fail later in `get_source` or
  `view_domain`, which then could not show the model at all.
- `create_domain` throws an `ArgumentError` for a cell width, size or resolution that is
  not positive and finite, or an axis without a cell. Zero, negative or NaN widths used
  to give nodes out of order, and with them silently wrong cells.
- `create_circular_loop_source!` and `create_rectangular_loop_source!` throw an
  `ArgumentError` on a grid with fewer than three nodes in a direction, where the range
  `nodes[1] … nodes[end-2]` is empty, instead of a `BoundsError`.
- `interpolate` of edge, facet and dual-node quantities, and with it the field view of
  `view_domain`, failed on a domain one cell thick along an axis ("Need at least two
  dual nodes"). The quantities centred in that cell are now constant along it.
- `get_position_of_index` threw a `BoundsError` whose message read "attempt to access
  35-codeunit String". It now names the valid index range of the entity, also for dead
  entities such as a u-edge leaving the last node plane in u.
- Objects, charges and interpolation points exactly on the far boundary of the domain
  were rejected, or interpolated as outside, when the last node fell below the requested
  extent by rounding (ten cells of 0.1 m end at 0.9999999999999999 m). Positions are now
  compared with the nodes up to 10⁻⁹ of the domain length.
- `plot_nodal_values(...; logscale = true)` failed as soon as the slice contained zero or
  negative values, which a potential routinely has. These values are now left out of the
  heatmap.
- `interpolate(domain, DualNode(), …)` returns `NaN` with a warning outside the domain,
  like every other quantity. It used to extrapolate silently from the outermost cell.
- `interpolate` throws a `DimensionMismatch` if the data does not hold one value per
  entity, `Np` for nodes and `3Np` for edges and facets. A vector of the other kind used
  to be read silently, e.g. the u-components of an edge vector as nodal values, or to
  fail with an index error deep inside.
- The docstrings of `get_index_entity`, `get_reluctivity` and the dual-facet averaging
  were separated from their functions by a blank line and therefore not attached to
  them; they now show up in the help.

### Documentation

- `plot_nodal_values` documents its keywords. The docstrings no longer contain
  Documenter `[…](@ref)` links, which showed up literally in the help.
- The notebooks with PEC walls build the projection with `get_boundary_matrix` instead of
  a copy of it of their own.
- The round solids (sphere, cylinder, torus, helix) share one routine that fills the
  cells, instead of four copies of it; the cells filled are the same.

- Docstrings for `create_domain`, `FITDomain`, `create_brick!`, `create_sphere!`,
  `interpolate`, `get_curl`, `get_gradient`, `get_divergence`, `get_permittivity`,
  `get_conductivity`, `get_vacuum_cfl_time` (which is exactly the stability limit on an
  equidistant grid, so a step slightly below it is needed), `get_position_of_index`,
  the boundary functions, and the exported types: every exported name now has one.
- `CITATION.cff`: the ORCIDs are fields of their own, no longer part of the
  affiliation, and the release date, not a valid date before, is left for the tag.
- The coil, eddy-current and SPFD dosimetry notebooks use `create_circular_loop_source!`
  with `current` and `get_source`. The eddy-current notebook names `create_brick!`
  where it still said `create_cube!`, and `Domain_and_History.ipynb` re-creates its
  sphere where it was, instead of on the copper shield.
- New notebook `examples/snippets/Domain_and_History.ipynb`.
- New script `examples/snippets/View_Domain.jl`: a model with every object type in
  `view_domain`. The examples environment includes GLMakie for it.
- New script `examples/snippets/View_Field.jl`: a dielectric sphere between two plates,
  its potential, E and D in the field view of `view_domain`.
- New script `examples/snippets/View_Magnetostatics.jl`: an iron C-core electromagnet
  with an air gap and a helical coil, B, H and A in the field view.
- README: figures of the field view, rendered from these two examples by
  `examples/snippets/make_readme_figures.jl` into `docs/images`.
- New notebook `examples/static/E-Statics-Dielectric-Sphere.ipynb`: a dielectric sphere
  in a uniform field on four grids via `change_resolution`, against the analytical
  solution.
- The examples find the package through a relative path in `examples/Project.toml`, so
  they run from any clone. The first cell of each notebook sets up the environment, on
  Julia 1.10 as well as on later versions.

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
