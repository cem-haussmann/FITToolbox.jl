# make_readme_figures.jl
# Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
# Date: 22/10/2026
#
# Renders the figures of the README from View_Field.jl and View_Magnetostatics.jl:
# the 3D scene of view_domain with a field on a cut, without the panel. Run it again
# when the viewer or the examples change:
#     julia examples/snippets/make_readme_figures.jl [output directory]
# The default output directory is docs/images. It sets the controls of the viewer
# through its internals (_view_figure and the blocks of the figure), as a user would
# in the window, so it may need adapting when the panel changes.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
VERSION < v"1.11" && Pkg.develop(path = joinpath(@__DIR__, "..", ".."); io = devnull)   # [sources] needs Julia ≥ 1.11
Pkg.instantiate(; io = devnull)
using GLMakie
using FITToolbox

# each example in a module of its own, so their variables do not meet; README_FIGURES
# tells them not to open a window
module Sphere
const README_FIGURES = true
include(joinpath(@__DIR__, "View_Field.jl"))
end
module Magnet
const README_FIGURES = true
include(joinpath(@__DIR__, "View_Magnetostatics.jl"))
end

const V = Base.get_extension(FITToolbox, :FITToolboxMakieExt)
out = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "docs", "images")
mkpath(out)

axis_menu(fig) = only(b for b in fig.content if b isa Menu && b.options[] == ["x", "y", "z"])
select!(m, option) = (m.i_selected[] = findfirst(==(option), m.options[]))

# The viewer with a field on the cut along y at `pos`, seen from `direction`; with
# `cut`, the objects on the camera's side of the plane are cut away. Returns the image
# of the 3D scene alone.
function field_scene(domain, fields; field, quantity, pos, log = false, arrows = true,
                     cut = false, hide = String[], direction = (0.45, 1.0, 0.5), distance = 1.9)
    fig, st = V._view_figure(domain, fields...; size = (1400, 750))
    f = st.field
    select!(axis_menu(fig), "y")
    slider = only(b for b in fig.content if b isa Slider && b.range[] == domain.nodes_v)
    set_close_to!(slider, pos)
    f.on.checked[] = true
    select!(f.field_menu, field)
    select!(f.quantity_menu, quantity)
    f.log[] = log
    f.arrows[] = arrows
    st.cut[] = cut
    for e in list_objects(domain).entries
        e.name in hide && (st.shown[e.obj_id][] = false)
    end
    # the camera: towards the middle of the domain, from `direction`, `distance` times
    # its diagonal away
    lo = (domain.nodes_u[1], domain.nodes_v[1], domain.nodes_w[1])
    hi = (domain.nodes_u[end], domain.nodes_v[end], domain.nodes_w[end])
    mid = (lo .+ hi) ./ 2
    dir = direction ./ sqrt(sum(abs2, direction))
    eye = mid .+ distance * sqrt(sum(abs2, hi .- lo)) .* dir
    cam = Makie.cameracontrols(st.ls.scene)
    update_cam!(st.ls.scene, cam, Vec3f(eye), Vec3f(mid), Vec3f(0, 0, 1))
    # the whole figure, cut down to the scene (saving the scene itself crops wrongly)
    img = colorbuffer(fig; px_per_unit = 1)
    (x0, y0), (w, h) = Tuple(st.ls.scene.viewport[].origin), Tuple(st.ls.scene.viewport[].widths)
    H = size(img, 1)
    return trim(img[H - (y0 + h) + 1:H - y0, x0 + 1:x0 + w])
end

# The image without its white margin, keeping `pad` pixels
function trim(img; pad = 10)
    used = [c != RGBf(1, 1, 1) for c in RGBf.(img)]
    rows, cols = findall(any(used; dims = 2)[:]), findall(any(used; dims = 1)[:])
    isempty(rows) && return img
    return img[max(1, first(rows) - pad):min(size(img, 1), last(rows) + pad),
               max(1, first(cols) - pad):min(size(img, 2), last(cols) + pad)]
end

# the dielectric sphere: |E| and its direction, the sphere hidden to show the field inside
img = field_scene(Sphere.domain, Sphere.fields; field = "E (V/m)", quantity = "magnitude",
                  pos = 0.5, hide = ["sphere"])
save(joinpath(out, "view_field_sphere.png"), img)

# the electromagnet: |B| on a log scale and its direction, in the middle of the core
img = field_scene(Magnet.domain, Magnet.fields; field = "B (T)", quantity = "magnitude",
                  pos = 0.1, log = true, cut = true)
save(joinpath(out, "view_field_electromagnet.png"), img)

println("written to ", abspath(out))
