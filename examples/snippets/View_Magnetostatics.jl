# View_Magnetostatics.jl
# Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
# Date: 22/10/2026
#
# Interactive 3D view of a magnetostatic field with view_domain: an electromagnet. An
# iron C-core (μ_r = 1000, legs 3 cm × 3 cm, 20 cm outside) with a 10 mm air gap in its
# right leg, and a coil of 6 turns carrying 10 A around its left leg.
#
# The magnetic circuit estimate for the gap, B ≈ μ0 N I / (g + l_Fe / μ_r) = 7.1 mT,
# neglects fringing and leakage; here the gap field is about 6 mT, and the left leg
# carries more flux than the right one: part of it closes across the window.
#
# Run from the repository root:
#     julia examples/snippets/View_Magnetostatics.jl
# The script waits until the window is closed. With `julia -i …` the REPL stays open
# instead. view_domain needs GLMakie (a window, as here) or WGLMakie (in a notebook).
#
# Things to try in the window:
# - tick "Field", and set the cut to y at 0.1 m: the plane through the middle of the
#   core; tick "cut" as well, to see the core and the coil cut open
# - B as "magnitude" with "log scale", and "arrows": the flux runs round the core,
#   fringes at the gap and leaks across the window
# - "show z": up the left leg, down the right one
# - H instead of B: small in the iron, large in the gap
# - "combined" view: the core with its window and gap cut out; "voxels": as the
#   solver sees it, and the coil on its grid edges

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
VERSION < v"1.11" && Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)   # [sources] needs Julia ≥ 1.11
Pkg.instantiate(; io = devnull)
using GLMakie
using FITToolbox
using LinearAlgebra
using IterativeSolvers

domain = create_domain([0.3, 0.2, 0.3], [0.005, 0.005, 0.005]; units = "m")

# the core: a brick, with the window and the gap cut out of it
create_brick!(domain, 0.05, 0.085, 0.05, 0.20, 0.03, 0.20; μ_r = 1000.0, name = "core", color = "gray")
create_brick!(domain, 0.08, 0.085, 0.08, 0.14, 0.03, 0.14; color = "background", name = "window")
create_brick!(domain, 0.22, 0.085, 0.145, 0.03, 0.03, 0.01; color = "background", name = "gap")

# the coil: a helix around the left leg (x 0.05 … 0.08 m), radius 3 cm, pitch 15 mm,
# starting and ending at -x, closed by a straight return wire at x = 15 mm
N, I_coil = 6, 10.0
c, ρ, pitch, z0 = (0.065, 0.1), 0.03, 0.015, 0.105
helix = [(c[1] + ρ * cos(θ), c[2] + ρ * sin(θ), z0 + pitch * (θ - π) / 2π)
         for θ in range(π, π + 2π * N; length = 32N + 1)]
lead = [(0.015, 0.1, helix[end][3]), (0.015, 0.1, helix[1][3])]
create_polyline_source!(domain, [helix; lead]; closed = true, current = I_coil, name = "coil")
j = get_source(domain, "coil")

# Magnetostatics Cᵀ M_ν C a = j. On the walls the tangential a is zero, so no flux
# leaves the box (B·n = 0); the natural condition of Cᵀ M_ν C, magnetic walls, would
# draw flux out of the core this close to it. CG solves the singular but consistent
# system, since the coil's j is free of divergence.
C  = get_curl(domain, Primal())
Mν = get_reluctivity(domain)
keep = ones(3 * domain.Np)
keep[get_all_tangential_boundary_indices(domain)] .= 0
P = Diagonal(keep)
K = P * (C' * Mν * C) * P + (I - P)
a, history = cg(K, P * j; reltol = 1e-8, maxiter = 20_000, log = true)
println("CG: converged = $(history.isconverged) after $(history.iters) iterations")

b = C * a            # facet fluxes:        PrimalFacet → B in T
h = Mν * b           # dual edge voltages:  DualEdge    → H in A/m
println("B in the gap: ", round(1e3 * interpolate(domain, PrimalFacet(), b, 0.235, 0.1, 0.15)[3];
                                 digits = 2), " mT")

list_objects(domain) |> display
fields = ("B (T)" => (b, PrimalFacet()), "H (A/m)" => (h, DualEdge()), "A (Vs/m)" => (a, PrimalEdge()))
# make_readme_figures.jl includes this file for the model and the fields alone
if !@isdefined(README_FIGURES)
    screen = display(view_domain(domain, fields...))
    # run as a script, Julia would exit here and close the window: wait until it is closed
    isinteractive() || wait(screen)
end
