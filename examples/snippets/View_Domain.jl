# View_Domain.jl
# Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
# Date: 06/10/2026
#
# Interactive 3D view of a model with view_domain.
#
# Run from the repository root:
#     julia examples/snippets/View_Domain.jl
# The script waits until the window is closed. With `julia -i …` the REPL stays open
# instead, so the domain can still be changed and viewed again.
#
# view_domain needs a Makie backend with a depth buffer: GLMakie (a window, as here)
# or WGLMakie (in a notebook). CairoMakie cannot draw overlapping 3D objects.
#
# Things to try in the window:
# - View "voxels": the cells each object occupies on the grid
# - Cut along y at 0.5 m: the core inside the ball, as modelled and as voxels,
#   and the loop between shield and ball
# - "grid in cut plane", and the loop "on grid edges" with "current direction"
# - click an object to hide it, "Show all" to bring it back

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
VERSION < v"1.11" && Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)   # [sources] needs Julia ≥ 1.11
Pkg.instantiate(; io = devnull)
using GLMakie
using FITToolbox

domain = create_domain([1.0, 1.0, 1.0], [0.05, 0.05, 0.05]; units = "m")

# stacked along z, which points up in the view: the shield as a horizontal plate at
# the bottom, the loop lying flat above it, the ball on top, so the view from above
# shows all of them; a rod stands on the shield beside them

# a plate from z = 0.1 m to 0.2 m
create_brick!(domain, 0.1, 0.1, 0.1, 0.8, 0.8, 0.1; units = "m", σ = 5.8e7,
              name = "shield", color = "copper")
# the ball from z = 0.4 m to 0.7 m
create_sphere!(domain, 0.5, 0.5, 0.55, 0.15; units = "m", ε_r = 4.0, name = "ball")
# a second sphere inside the first: hidden unless you cut or hide "ball"
create_sphere!(domain, 0.5, 0.5, 0.55, 0.08; units = "m", ε_r = 9.0,
               name = "core", color = "red")
# the loop in the plane z = 0.3 m, between the shield and the ball
create_circular_loop_source!(domain, 0.1, 0.5, 0.5, 0.3, DirZ(); units = "m",
                             current = 1000.0, name = "coil")
# an aluminium rod standing on the shield next to the ball, z = 0.2 m to 0.8 m:
# the base face is at (0.2, 0.75, 0.2), the rod extends 0.6 m along +z
create_cylinder!(domain, 0.2, 0.75, 0.2, 0.1, 0.6, DirZ(); units = "m", σ = 3.5e7,
                 name = "rod", color = "aluminium")

list_objects(domain) |> display
screen = display(view_domain(domain))
# run as a script, Julia would exit here and close the window: wait until it is closed
isinteractive() || wait(screen)
