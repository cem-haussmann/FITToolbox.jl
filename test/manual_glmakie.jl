# manual_glmakie.jl
# Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
# Date: 30/10/2026
#
# view_domain drawn by GLMakie: what the test suite cannot check, since CI has no GPU and
# CairoMakie builds the figure but cannot draw it. Not part of runtests.jl; run it by
# hand after changes to the viewer, on a machine with a display:
#
#     julia -t 2 test/manual_glmakie.jl [output folder]
#
# The windows stay hidden. It checks that every view is drawn and differs as it should,
# that hiding objects changes the picture and Show all restores it, that a click picks
# the object under the mouse (one plot per object, and merged plots), that "play" steps
# through a time series at the speed chosen, and that the colour limits reach the picture; it prints the
# time per frame with 1000 objects and saves a PNG of each view for a look by eye. It
# uses the environment of the examples, which has GLMakie.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "examples"); io = devnull)
VERSION < v"1.11" && Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)   # [sources] needs Julia ≥ 1.11
Pkg.instantiate(; io = devnull)
using FITToolbox, GLMakie, Test, Printf

GLMakie.activate!(visible = false)
const V = Base.get_extension(FITToolbox, :FITToolboxMakieExt)
const OUT = isempty(ARGS) ? mkpath(joinpath(tempdir(), "fittoolbox_glmakie")) : mkpath(ARGS[1])
quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

# --- helpers --------------------------------------------------------------------

# the picture, and saved under `name`: a copy, as colorbuffer reuses its buffer
function frame(fig, name = nothing)
    img = copy(Makie.colorbuffer(fig))
    name === nothing || save(joinpath(OUT, name * ".png"), img)
    return img
end
# the pixels in which two pictures differ, and those not near white (something drawn)
ndiff(a, b) = count(i -> a[i] != b[i], eachindex(a, b))
ninked(img) = count(c -> min(Makie.Colors.red(c), Makie.Colors.green(c), Makie.Colors.blue(c)) < 0.9, img)

# what a click at pixel `px` of the scene picks: the obj_id, or nothing
function picked(st, px)
    plt, idx = pick(st.ls.scene, Point2f(px))
    return V._picked_id(st.plot_ids, plt, idx)
end
# pixels on a lattice over the scene, for picking
function lattice(st; n = 25)
    vp = st.ls.scene.viewport[]
    o, w = minimum(vp), widths(vp)
    return [Point2f(o[1] + w[1] * i / (n + 1), o[2] + w[2] * j / (n + 1)) for i in 1:n, j in 1:n]
end

show_all!(fig) = only(b for b in fig.content if b isa Button && b.label[] == "Show all").clicks[] = 1

# a model with every kind of object, few enough for a plot per object
function every_kind()
    d = create_domain([1.0, 1.0, 1.0], [0.02, 0.02, 0.02])
    quiet() do
        create_brick!(d, 0.1, 0.1, 0.05, 0.8, 0.8, 0.04; σ = 5.8e7, name = "plate", color = "copper")
        create_sphere!(d, 0.5, 0.5, 0.45, 0.2; ε_r = 4.0, name = "ball")
        create_sphere!(d, 0.5, 0.5, 0.45, 0.08; color = "background", name = "hole")
        create_cylinder!(d, 0.2, 0.8, 0.1, 0.04, 0.6, DirZ(); μ_r = 100.0, name = "rod")
        create_torus!(d, 0.5, 0.5, 0.8, 0.25, 0.03, DirZ(); σ = 1e6, name = "ring")
        create_helix!(d, 0.8, 0.2, 0.1, 0.08, 0.015, 0.06, 5, DirZ(); σ = 5.8e7, name = "coil")
        create_circular_loop_source!(d, 0.5, 0.5, 0.7, 0.3, DirZ(); name = "loop")
        create_rectangular_loop_source!(d, 0.1, 0.6, 0.3, 0.3, 0.3, DirX(); name = "frame")
        create_polyline_source!(d, [(0.1, 0.1, 0.9), (0.9, 0.2, 0.9), (0.9, 0.9, 0.6)]; name = "wire")
        create_electric_charge!(d, 0.3, 0.3, 0.9; charge = 1.0, name = "q")
    end
    return d
end

# many small spheres in three colours, and a few other objects: merged plots
function many(n)
    d = create_domain([1.0, 1.0, 1.0], [0.01, 0.01, 0.01])
    quiet() do
        for k in 1:n
            create_sphere!(d, 0.05 + 0.9 * ((k-1) % 30) / 29, 0.05 + 0.9 * (((k-1) ÷ 30) % 30) / 29,
                           0.1 + 0.4 * ((k-1) ÷ 900), 0.012; σ = 1.0,
                           color = ("copper", "aluminium", "gold")[mod1(k, 3)])
        end
        create_circular_loop_source!(d, 0.5, 0.5, 0.75, 0.2, DirZ(); name = "loop")
        create_brick!(d, 0.3, 0.3, 0.6, 0.4, 0.4, 0.05; color = "background", name = "window")
        create_electric_charge!(d, 0.5, 0.5, 0.9; charge = 1.0, name = "q")
    end
    return d
end

# --- the checks -------------------------------------------------------------------

@testset "view_domain with GLMakie" begin

@testset "view_domain accepts GLMakie" begin
    @test quiet(() -> view_domain(every_kind())) isa Figure
end

@testset "every view, one plot per object" begin
    d = every_kind()
    fig, st = quiet(() -> V._view_figure(d))
    @test st.merged === nothing
    defined = frame(fig, "1_defined")
    @test ninked(defined) > 10_000
    st.combined_mode[] = true
    combined = frame(fig, "2_combined")
    @test ndiff(combined, defined) > 1000
    st.combined_mode[] = false; st.voxel_mode[] = true
    voxels = frame(fig, "3_voxels")
    @test ndiff(voxels, defined) > 1000 && ndiff(voxels, combined) > 1000
    st.cut[] = true
    cut = frame(fig, "4_voxels_cut")
    @test ndiff(cut, voxels) > 1000
    st.cut[] = false; st.voxel_mode[] = false
    # hidden: another picture; Show all: the same as before
    st.shown[2][] = false
    @test ndiff(frame(fig), defined) > 1000
    show_all!(fig)
    @test ndiff(frame(fig), defined) == 0
end

@testset "a click picks the object under the mouse" begin
    for (label, d) in (("a plot per object", every_kind()), ("merged plots", many(200)))
        fig, st = quiet(() -> V._view_figure(d))
        frame(fig)
        ids = [picked(st, px) for px in lattice(st)]
        hits = filter(!isnothing, ids)
        # objects under some of the pixels, each one that exists
        @test length(unique(hits)) >= 3
        @test all(id -> haskey(st.shown, id), hits)
        # hide one that was hit: where it was picked, it is not any more
        id = first(hits)
        at = [px for (px, i) in zip(lattice(st), ids) if i == id]
        st.shown[id][] = false
        frame(fig)
        @test all(px -> picked(st, px) != id, at)
        println("  picking, $label: $(length(unique(hits))) objects hit on a 25 × 25 lattice")
    end
end

@testset "1000 objects: merged, and the time per frame" begin
    d = many(1000)
    quiet(() -> V._view_figure(many(30))) |> first |> frame        # compile
    fig, st = quiet(() -> V._view_figure(d))
    @test length(st.ls.scene.plots) < 40
    defined = frame(fig, "5_many_defined")
    cam = Makie.cameracontrols(st.ls.scene)
    t = @elapsed for k in 1:20
        θ = 2π * k / 20
        Makie.update_cam!(st.ls.scene, cam, Vec3f(3cos(θ), 3sin(θ), 2), Vec3f(0.5), Vec3f(0, 0, 1))
        Makie.colorbuffer(fig)
    end
    @printf("  1003 objects in %d plots: %.1f ms per frame\n", length(st.ls.scene.plots), 1000t / 20)
    V._reset_view!(st.ls.scene)
    for k in 1:3:300; st.shown[k][] = false; end
    @test ndiff(frame(fig, "6_many_hidden"), defined) > 1000
    show_all!(fig)
    st.voxel_mode[] = true
    @test ndiff(frame(fig, "7_many_voxels"), defined) > 1000
end

@testset "a field: colour limits and play" begin
    d = create_domain([1.0, 1.0, 1.0], [0.02, 0.02, 0.02])
    create_sphere!(d, 0.5, 0.5, 0.5, 0.2; ε_r = 4.0, name = "ball")
    Φ = [10x - 5 + 2exp(-((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2) / 0.01)
         for z in d.nodes_w for y in d.nodes_v for x in d.nodes_u]
    snaps = [k .* Φ for k in 1:5]
    fig, st = V._view_figure(d, "Φ (V)" => (Φ, PrimalNode()),
                             "Φ(t) (V)" => (snaps, PrimalNode(), collect(1e-9:1e-9:5e-9)))
    f = st.field
    f.on.checked[] = true
    st.cut[] = true
    auto = frame(fig, "8_field")
    f.min_box.stored_string[] = "0"
    set_close_to!(f.range_slider, 0.0, 3.0)
    @test all(isapprox.(f.colorrange[], (0.0, 3.0); rtol = 1e-6))
    @test ndiff(frame(fig, "9_field_limits"), auto) > 1000
    # play: in a window (hidden, but open) the steps advance until paused
    f.field_menu.i_selected[] = 2
    k0 = f.step[]
    f.play.clicks[] = 1
    t0 = time()
    while f.step[] == k0 && time() - t0 < 5
        sleep(0.05)
    end
    @test f.step[] != k0 && f.playing[]
    f.play.clicks[] = 1
    sleep(0.3)
    @test !f.playing[]
    # the speed chosen: steps in 2 s at 1 and at 10 snapshots per second (a step takes
    # its time, so 10 / s may fall short of 20 on a slow machine); pause stops at once
    steps = Dict{Float64, Int}()
    for (i, s) in ((2, 1.0), (5, 10.0))
        f.speed_menu.i_selected[] = i
        n = Ref(0)
        counter = on(_ -> n[] += 1, f.step)
        f.play.clicks[] = 1
        sleep(2.0)
        f.play.clicks[] = 1
        off(counter)
        steps[s] = n[]
    end
    println("  play: $(steps[1.0]) steps in 2 s at 1 / s, $(steps[10.0]) at 10 / s")
    @test 2 <= steps[1.0] <= 3
    @test steps[10.0] >= 10
    sleep(0.2)
    @test !f.playing[]
end

end

println("pictures in ", OUT)
