# View_Field.jl
# Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
# Date: 22/10/2026
#
# Interactive 3D view of a model and its field with view_domain: a dielectric sphere
# (ε_r = 4, radius 15 cm) between two plates 1 m apart at 0 V and 1 V, i.e. in a
# uniform field E0 = 1 V/m along -z. Inside the sphere the field is uniform as well,
# 3 / (ε_r + 2) E0 = 0.5 V/m; outside it bends around the sphere.
#
# Run from the repository root:
#     julia examples/snippets/View_Field.jl
# The script waits until the window is closed. With `julia -i …` the REPL stays open
# instead, so the fields can be changed and viewed again.
#
# view_domain needs a Makie backend with a depth buffer: GLMakie (a window, as here)
# or WGLMakie (in a notebook). CairoMakie cannot draw overlapping 3D objects.
#
# Things to try in the window:
# - tick "Field", and set the cut to y at 0.5 m: the plane runs through the sphere
#   along the field; click the sphere to hide it
# - E with "show z": -0.5 V/m inside the sphere, stronger above and below it, weaker
#   at its sides; "fit colours to plane" and "magnitude" use the colour map better
# - "arrows": the field lines bending around the sphere; "scaled length" adds the
#   strength
# - D instead of E: its normal component is continuous across the sphere's surface,
#   E's is not
# - Φ: the potential, 0 V at the bottom plate and 1 V at the top
# - move the slider of the cut: the plane follows it

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
VERSION < v"1.11" && Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)   # [sources] needs Julia ≥ 1.11
Pkg.instantiate(; io = devnull)
using GLMakie
using FITToolbox
using LinearAlgebra

domain = create_domain([1.0, 1.0, 1.0], [0.02, 0.02, 0.02]; units = "m")
create_sphere!(domain, 0.5, 0.5, 0.5, 0.15; units = "m", ε_r = 4.0, name = "sphere")

# Laplace's equation L Φ = 0 with L = Gᵀ M_ε G. The plates z = 0 and z = 1 m are
# Dirichlet boundaries at 0 V and 1 V; the other four faces are Neumann (no flux).
# R keeps the free nodes, D the fixed ones; Φ0 holds the plate voltages.
G  = get_gradient(domain, Primal())
Mε = get_permittivity(domain)
L  = G' * Mε * G
plates = [(DirZ(), Negative()), (DirZ(), Positive())]
R  = get_boundary_matrix(domain, plates, NodalComponent())
D  = Diagonal(ones(domain.Np)) - R
Φ0 = zeros(domain.Np)
Φ0[get_custom_boundary_indices(domain, [(DirZ(), Positive())], NodalComponent())] .= 1.0
Φ  = cholesky(D + R * L * R; check = true) \ (-(R * (L * Φ0))) + Φ0

# the field quantities on their entities, as the view expects them
e = -G * Φ           # edge voltages:      PrimalEdge → E in V/m
d = Mε * e           # dual facet fluxes:  DualFacet  → D in C/m²

list_objects(domain) |> display
fields = ("Φ (V)" => (Φ, PrimalNode()), "E (V/m)" => (e, PrimalEdge()), "D (C/m²)" => (d, DualFacet()))
# make_readme_figures.jl includes this file for the model and the fields alone
if !@isdefined(README_FIGURES)
    screen = display(view_domain(domain, fields...))
    # run as a script, Julia would exit here and close the window: wait until it is closed
    isinteractive() || wait(screen)
end
