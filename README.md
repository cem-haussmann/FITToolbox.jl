# FITToolbox.jl

[![DOI](https://zenodo.org/badge/1355999732.svg)](https://doi.org/10.5281/zenodo.22282722)

A Julia implementation of the **Finite Integration Technique** (FIT) for
computational electromagnetics on structured, optionally non-equidistant grids.

The package provides the discrete topological operators, the material matrices,
boundary conditions, and field reconstruction on the staggered primal–dual grid
pair — the building blocks from which electrostatic, magnetostatic and (in time)
full-wave solvers are assembled.

> **Status:** early development. The API is not yet stable and features are being
> added as they are needed. Feedback and issues are welcome.

## Why FIT

The Finite Integration Technique was introduced by Thomas Weiland in 1977 and
discretises Maxwell's equations in their **integral** form on a pair of staggered
grids. Field quantities are not sampled at points but integrated along edges and
over facets, which gives the method several properties that are convenient in
practice:

- **The topology is exact.** The discrete curl, gradient and divergence contain
  only ±1 entries and satisfy `curl grad = 0` and `div curl = 0` *identically*, to
  machine precision, on any grid. Charge and energy conservation follow from the
  operators themselves rather than from careful discretisation.
- **All approximation lives in one place.** Geometry and material enter only
  through the diagonal material matrices `M_ε`, `M_ν`, `M_σ`. Everything else is
  metric-free. This makes the method easy to reason about and easy to extend —
  anisotropy, dispersion or a graded mesh change the material matrices and nothing
  else.
- **The equations translate directly.** The discrete Poisson problem is
  `Gᵀ M_ε G Φ = q`; magnetostatics is `Cᵀ M_ν C â = ĵ`. The matrix form mirrors the
  continuous operators one-to-one, so the step from theory to code is short.
- **It generalises FDTD.** On a uniform Cartesian grid with explicit time stepping,
  FIT and the Yee scheme are computationally equivalent — but the matrix notation
  makes stability, consistency and conservation properties directly analysable.

## Installation

```julia
julia> ]
pkg> add https://github.com/cem-haussmann/FITToolbox.jl.git
```

Requires Julia 1.10 (the current long-term-support release) or later.

The test suite covers the operator identities, the material matrix values on both
uniform and graded grids, the boundary projections and the interpolation:

```julia
pkg> test FITToolbox
```

## Example: the potential of a point charge

```julia
using FITToolbox, LinearAlgebra

# 100 m cube at 2 m resolution, vacuum
domain = create_domain([100.0, 100.0, 100.0], [2.0, 2.0, 2.0];
                       units = "m", σ = 0, ε_r = 1, μ_r = 1)

# A proton at the centre. Charge in FIT is already volume-integrated, so the
# value is written directly into the dual volume that contains it.
q = zeros(domain.Np)
q[get_index_entity(domain, DualVolume(), 50.0, 50.0, 50.0)] = 1.602e-19

# Discrete Poisson operator. The ghost matrix removes the "dead" edges in the
# last plane of each direction, which have no geometric counterpart.
G = get_ghost_matrix(domain) * get_gradient(domain, Primal())
L = G' * get_permittivity(domain) * G

# Ground the outer boundary: R keeps the interior, D supplies the diagonal
# entries for the rows it zeroes.
faces = [(n, s) for n in (X(), Y(), Z()) for s in (Positive(), Negative())]
R = get_boundary_matrix(domain, faces, NodalComponent())
D = Diagonal(ones(domain.Np)) - R

Φ = cholesky(D + R*L*R; check = true) \ (R * q)

# Sample the potential anywhere, not just on nodes
interpolate(domain, PrimalNode(), Φ, 60.0, 50.0, 50.0; units = "m")
```

Non-equidistant grids are built by passing the cell widths directly, which is
useful for refining a region of interest or for pushing the boundary far away at
modest cost:

```julia
edges = vcat(20.0, 10.0, 5.0, fill(2.0, 50), 5.0, 10.0, 20.0)
domain = create_domain(edges, edges, edges; units = "m")
```

## Examples

Worked notebooks are in `examples/`, ordered roughly by increasing difficulty.
Each explains its own physics and numerics rather than presenting a finished
recipe, and each compares against an analytical solution where one exists.

**Electrostatics** (`examples/static/`)

| Notebook | Topic |
|---|---|
| `E-Statics-PC-Neumann-Dirichlet.ipynb` | Point charge with homogeneous Neumann and Dirichlet conditions. Why the Neumann system is both singular and inconsistent, how gauge fixing and charge neutralisation differ, and what each boundary does to the field. |
| `E-Statics-PC-Robin-Graded.ipynb` | Robin conditions, which prescribe the far-field decay instead of the value or the flux, and domain extension by graded cells. Compares the cost of each against the accuracy it buys. |

**Magnetostatics and eddy currents** (`examples/static/`, `examples/lf/`)

| Notebook | Topic |
|---|---|
| `static/M-Statics-Coil.ipynb` | Curl–curl system for a circular coil. The gauge freedom spans roughly `N_p` dimensions, so no finite set of pinned nodes removes it; an iterative solver is used instead. Contrasts electric and magnetic wall conditions. |
| `lf/M-Quasi-Statics-Coil-Shielding.ipynb` | Time-harmonic eddy currents in a conducting plate at 50 Hz. The complex symmetric system, why skin depth rather than conductivity decides whether shielding is visible. |


**Time domain** (`examples/hf/`)

| Notebook | Topic |
|---|---|
| `FDTD-Dipole-PEC.ipynb` | Explicit leapfrog time stepping for a radiating dipole in a closed metal box. The Courant limit, the staggering in time, and the reflections a PEC wall produces. |
| `FDTD-Dipole-Sphere-PEC.ipynb` | Scattering from a conducting sphere. Exponential integration of the conductive term, which stays stable at the vacuum time step where the usual semi-implicit update fails by seven orders of magnitude. |
| `FDTD-Dipole-ABC.ipynb` | A matched absorbing layer: graded electric conductivity together with a magnetic loss chosen to preserve the wave impedance, so the boundary attenuates instead of reflecting. |

**Dosimetry** (`examples/lf/`)

| Notebook | Topic |
|---|---|
| `Dosimetry-SPFD.ipynb` | Induced field in a nine-layer conducting sphere, computed with the scalar-potential finite-difference scheme. Why the problem splits into a curl–curl solve for the source and a Poisson solve for the charge, and what the split costs on a voxel mesh. Reproduces the benchmark of Conchin Gubernati et al. (2022). |

This last notebook is worth singling out. It is a reproduction of a published
benchmark, and its result is a negative one: the source field alone matches the
analytical reference to 1–3 %, while adding the physically required charge
correction degrades the per-layer maxima to 12–33 %. The benchmark is constructed so
that the exact correction vanishes by symmetry, so everything the correction adds is
staircasing — which is what the source paper set out to demonstrate. The notebook
therefore validates the implementation and illustrates the limitation of Cartesian
meshes for curved tissue boundaries at the same time.

The time-domain solver is verified against the closed-form field of a Hertzian
dipole, which is exact at every radius rather than only in the far zone. The
computed envelope agrees to within 1–3 % across the radiation region on both a
uniform and a graded grid; the departures at very small radius (where a one-cell
source is not a point) and near the absorbing layer are expected and explained in
the notebook.

## Intended audience

The package is written for two groups. For **students**, the aim is that the code
should read like the theory: the operators are named after what they are, the
matrix expressions mirror the equations in the lecture notes, and the example
notebooks explain each step rather than presenting a finished recipe. For
**academic research**, the aim is a clean, hackable basis for method development —
new boundary conditions, material models or solver strategies can be tried without
fighting a framework.

It is not a production solver. For large industrial models, established codes are
faster and far better validated.

## On the use of large language models

The original code predates the availability of capable LLMs and was written by
hand. For the version presented here, LLMs were used for restructuring and
tidying the code, for improving readability and consistency, for a large part of
the comments in both the source and the example notebooks, and for this README.

All physics, all numerical methods and all design decisions remain the author's.
Every LLM-suggested change was reviewed before being adopted. Users should nonetheless
treat the code as they would any research software: verify results against known
solutions before trusting them.

## References

The method:

1. T. Weiland, "A discretization method for the solution of Maxwell's equations
   for six-component fields", *Electronics and Communications (AEÜ)*, vol. 31,
   no. 3, pp. 116–120, 1977. — the original paper introducing FIT.
2. T. Weiland, "On the unique numerical solution of Maxwellian eigenvalue problems
   in three dimensions", *Particle Accelerators*, vol. 17, pp. 227–242, 1985.
3. T. Weiland, "Time domain electromagnetic field computation with finite
   difference methods", *International Journal of Numerical Modelling*, vol. 9,
   pp. 295–319, 1996.
4. M. Clemens and T. Weiland, "Discrete electromagnetism with the finite
   integration technique", *Progress In Electromagnetics Research (PIER)*,
   vol. 32, pp. 65–87, 2001. doi:10.2528/PIER00080103 — the standard modern
   review, and the best starting point for the matrix notation used here.
5. R. Schuhmann and T. Weiland, "Conservation of discrete energy and related laws
   in the finite integration technique", *PIER*, vol. 32, pp. 301–316, 2001.

For comparison with the closely related Yee scheme:

6. K. S. Yee, "Numerical solution of initial boundary value problems involving
   Maxwell's equations in isotropic media", *IEEE Transactions on Antennas and
   Propagation*, vol. 14, no. 3, pp. 302–307, 1966.

## Citing

If this package is useful in your work, a citation is appreciated. See
`CITATION.cff`, or use the "Cite this repository" button in the GitHub sidebar.

## License

MIT — see `LICENSE`. You are free to use, modify and redistribute the code,
including commercially, provided the copyright notice is retained.

## Author

Norman Haußmann — Chair of Electromagnetic Theory, University of Wuppertal
<haussmann@uni-wuppertal.de>
Markus Clemens — Chair of Electromagnetic Theory, University of Wuppertal
<clemens@uni-wuppertal.de>