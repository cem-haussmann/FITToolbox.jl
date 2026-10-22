#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 01/10/2026
#
# Tests for object creation (CreateBrick.jl, CreateSphere.jl). Included from runtests.jl.

using FITToolbox, Test

@testset "create_brick!" begin
    d = create_domain([6.0, 8.0, 10.0], [2.0, 2.0, 2.0])
    create_brick!(d, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0; ε_r=5.0)
    @test all(d.material[1:2, 1:2, 1:2, 2] .== 5.0)   # cells inside the brick
    @test d.material[3, 1, 1, 2] == 1.0               # outside keeps background
end

@testset "create_sphere!" begin
    d = create_domain([6.0, 8.0, 10.0], [2.0, 2.0, 2.0])
    create_sphere!(d, 3.0, 4.0, 5.0, 2.0; ε_r=5.0, μ_r=2.0)
    @test any(d.material[:,:,:,2] .== 5.0)     # something was filled
    @test any(d.material[:,:,:,2] .== 1.0)     # but not everything
end

@testset "create_cylinder!" begin
    # 2 cm cells; a rod of radius 10 cm standing on z = 20 cm, 40 cm long
    mk() = create_domain([1.0, 1.0, 1.0], [0.02, 0.02, 0.02])
    d = mk()
    id = create_cylinder!(d, 50, 50, 20, 10, 40, DirZ(); units = "cm", ε_r = 4.0,
                          name = "rod", color = "copper")
    @test id == 1
    filled = d.material[:, :, :, 2] .== 4.0

    # along the axis exact: the end faces are on nodes, z = 0.2 … 0.6 m = layers 11:30
    layers = [any(filled[:, :, k]) for k in 1:size(filled, 3)]
    @test findall(layers) == collect(11:30)
    # across it: the same staircase circle in every layer
    @test all(k -> filled[:, :, k] == filled[:, :, 11], 11:30)
    # volume close to π R² H
    @test count(filled) * 0.02^3 ≈ π * 0.1^2 * 0.4 rtol = 0.02

    # recorded, and shown in the object table
    step = only(d._history.steps)
    @test step.object isa FITToolbox.Cylinder && step.name == "rod"
    @test occursin("base=(0.5, 0.5, 0.2) radius=0.1 height=0.4 axis=DirZ",
                   repr(MIME"text/plain"(), list_objects(d)))

    # the three axes give the same cylinder, turned: x ↔ z swaps the array dimensions
    x = mk(); create_cylinder!(x, 0.2, 0.5, 0.5, 0.1, 0.4, DirX(); ε_r = 4.0)
    z = mk(); create_cylinder!(z, 0.5, 0.5, 0.2, 0.1, 0.4, DirZ(); ε_r = 4.0)
    y = mk(); create_cylinder!(y, 0.5, 0.2, 0.5, 0.1, 0.4, DirY(); ε_r = 4.0)
    @test permutedims(x.material[:, :, :, 2], (3, 2, 1)) == z.material[:, :, :, 2]
    @test permutedims(y.material[:, :, :, 2], (1, 3, 2)) == z.material[:, :, :, 2]

    # the same 50 % rule as the sphere: a cylinder of height 2R and a sphere of radius
    # R fill the same cells in the middle layer
    s = mk(); create_sphere!(s, 0.5, 0.5, 0.5, 0.1; ε_r = 4.0)
    c = mk(); create_cylinder!(c, 0.5, 0.5, 0.4, 0.1, 0.2, DirZ(); ε_r = 4.0)
    @test s.material[:, :, 25, 2] == c.material[:, :, 25, 2]

    # invalid requests: nothing recorded
    e = mk()
    @test_throws "radius must be positive" create_cylinder!(e, 0.5, 0.5, 0.2, 0.0, 0.4, DirZ())
    @test_throws "height must be positive" create_cylinder!(e, 0.5, 0.5, 0.2, 0.1, -0.4, DirZ())
    @test_throws "outside the domain"      create_cylinder!(e, 0.5, 0.5, 1.2, 0.1, 0.4, DirZ())
    @test isempty(e._history.steps)

    # reaching past the boundary: clipped, with a warning
    @test_logs (:warn, r"clipped") create_cylinder!(e, 0.5, 0.5, 0.8, 0.1, 0.4, DirZ())
end

@testset "objects reaching the boundary exactly" begin
    # ten cells of 0.1 m end at 0.9999999999999999 m: objects up to 1.0 m must still fit
    d = create_domain([1.0, 1.0, 1.0], [0.1, 0.1, 0.1])
    @test d.nodes_u[end] < 1.0                                         # the rounding this is about
    @test create_brick!(d, 0.0, 0.0, 0.0, 1.0, 1.0, 0.1; ε_r = 2.0) == 1
    @test all(d.material[:, :, 1, 2] .== 2.0)                         # every cell of the bottom layer
    @test_logs create_cylinder!(d, 0.5, 0.5, 0.0, 0.2, 1.0, DirZ())   # up to the top: not clipped
    @test_logs (:warn, r"clipped") create_sphere!(d, 1.0, 0.5, 0.5, 0.2)   # centre on the boundary
    @test create_electric_charge!(d, 1.0, 1.0, 1.0; charge = 1e-9) == 4
    @test get_source(d, 4)[end] == 1e-9                               # the last node

    # still an error beyond it
    @test_throws "beyond the domain" create_brick!(d, 0.0, 0.0, 0.0, 1.001, 1.0, 0.1)
    @test_throws "outside the domain" create_electric_charge!(d, 1.001, 0.5, 0.5; charge = 1e-9)
end

@testset "create_torus!" begin
    # 2 cm cells; a ring lying flat at z = 50 cm, ring radius 25 cm, tube radius 8 cm,
    # i.e. 4 cells: thinner tubes make the staircase error several percent
    mk() = create_domain([1.0, 1.0, 1.0], [0.02, 0.02, 0.02])
    d = mk()
    id = create_torus!(d, 50, 50, 50, 25, 8, DirZ(); units = "cm", μ_r = 100.0,
                       name = "core", color = "darkgray")
    @test id == 1
    filled = d.material[:, :, :, 3] .== 100.0
    # volume close to 2π² R r²
    @test count(filled) * 0.02^3 ≈ 2π^2 * 0.25 * 0.08^2 rtol = 0.03
    # the hole is empty, the tube is filled, nothing above r from its plane
    @test !any(filled[24:27, 24:27, :])                            # around the axis
    @test filled[13, 25, 25] && filled[38, 25, 25]                 # x = 0.25 m and 0.75 m: in the tube
    @test !any(filled[:, :, 1:20]) && !any(filled[:, :, 31:end])   # |z - 0.5| > r
    # symmetric about its plane z = 0.5 m, the node between cell layers 25 and 26
    @test filled[:, :, 21:25] == filled[:, :, 30:-1:26]

    # recorded, and shown in the object table
    step = only(d._history.steps)
    @test step.object isa FITToolbox.Torus && step.name == "core"
    @test occursin("center=(0.5, 0.5, 0.5) R=0.25 r=0.08 axis=DirZ",
                   repr(MIME"text/plain"(), list_objects(d)))
    undo!(d)
    @test all(d.material[:, :, :, 3] .== 1.0)

    # the three axes give the same torus, turned
    x = mk(); create_torus!(x, 0.5, 0.5, 0.5, 0.2, 0.06, DirX(); ε_r = 4.0)
    y = mk(); create_torus!(y, 0.5, 0.5, 0.5, 0.2, 0.06, DirY(); ε_r = 4.0)
    z = mk(); create_torus!(z, 0.5, 0.5, 0.5, 0.2, 0.06, DirZ(); ε_r = 4.0)
    @test permutedims(x.material[:, :, :, 2], (3, 2, 1)) == z.material[:, :, :, 2]
    @test permutedims(y.material[:, :, :, 2], (1, 3, 2)) == z.material[:, :, :, 2]

    # invalid requests: nothing recorded
    e = mk()
    @test_throws "ring radius R must be positive" create_torus!(e, 0.5, 0.5, 0.5, 0.0, 0.05, DirZ())
    @test_throws "tube radius r must be positive" create_torus!(e, 0.5, 0.5, 0.5, 0.2, -0.05, DirZ())
    @test_throws "closes the hole"                create_torus!(e, 0.5, 0.5, 0.5, 0.2, 0.2, DirZ())
    @test_throws "outside the domain"             create_torus!(e, 0.5, 1.2, 0.5, 0.2, 0.05, DirZ())
    @test isempty(e._history.steps)

    # reaching past the boundary: clipped, with a warning; along the axis only r counts
    @test_logs (:warn, r"clipped") create_torus!(e, 0.5, 0.25, 0.5, 0.2, 0.06, DirZ())
    @test_logs create_torus!(e, 0.5, 0.5, 0.06, 0.2, 0.06, DirZ())   # touches z = 0: not clipped
    @test_logs (:warn, r"no cells were filled") create_torus!(e, 0.5, 0.5, 0.5, 0.2, 0.005, DirZ())
end

@testset "create_helix!" begin
    # 2 mm cells; a winding around z: radius 40 mm, wire radius 8 mm (4 cells), 3 turns
    # of 30 mm, starting at +x of the axis point (0.1, 0.1, 0.05)
    mk() = create_domain([0.2, 0.2, 0.2], [0.002, 0.002, 0.002])
    cell(d, p) = ntuple(i -> searchsortedlast((d.nodes_u, d.nodes_v, d.nodes_w)[i], p[i]), 3)
    d = mk()
    id = create_helix!(d, 100, 100, 50, 40, 8, 30, 3, DirZ(); units = "mm", σ = 5.8e7,
                       name = "coil", color = "copper")
    @test id == 1
    filled = d.material[:, :, :, 1] .== 5.8e7
    # volume close to the tube, π r² × length of the centre line, plus one sphere for the two caps
    L = 3 * hypot(2π * 0.04, 0.03)
    @test count(filled) * 0.002^3 ≈ π * 0.008^2 * L + 4π / 3 * 0.008^3 rtol = 0.03

    # the ends, at the centres of the round caps
    ends = get_helix_ends(d, "coil")
    @test all(ends.start .≈ (0.14, 0.1, 0.05)) && all(ends.stop .≈ (0.14, 0.1, 0.14))
    @test get_helix_ends(d, 1) == ends
    @test filled[cell(d, ends.start)...] && filled[cell(d, ends.stop)...]
    @test filled[cell(d, ends.start .- (0, 0, 0.006))...]       # the cap below the start
    # right-handed about +z: a quarter turn later at +y, a quarter pitch higher
    @test filled[cell(d, (0.1, 0.141, 0.0575))...] && !filled[cell(d, (0.1, 0.059, 0.0575))...]
    # empty on the axis and between the turns
    @test !any(filled[48:53, 48:53, :])
    @test !filled[cell(d, (0.141, 0.1, 0.065))...]

    # recorded, and shown in the object table
    step = only(d._history.steps)
    @test step.object isa FITToolbox.Helix && step.name == "coil"
    @test occursin("base=(0.1, 0.1, 0.05) R=0.04 r=0.008 pitch=0.03 turns=3.0 axis=DirZ start=0.0° right-handed",
                   repr(MIME"text/plain"(), list_objects(d)))

    # left-handed: the quarter turn goes to -y; start_angle turns the start
    l = mk(); create_helix!(l, 0.1, 0.1, 0.05, 0.04, 0.008, 0.03, 3, DirZ(); left_handed = true, ε_r = 4.0)
    @test l.material[cell(l, (0.1, 0.059, 0.0575))..., 2] == 4.0
    @test l.material[cell(l, (0.1, 0.141, 0.0575))..., 2] == 1.0
    a = mk(); id_a = create_helix!(a, 0.1, 0.1, 0.05, 0.04, 0.008, 0.03, 2.5, DirZ(); start_angle = 90)
    @test all(get_helix_ends(a, id_a).start .≈ (0.1, 0.14, 0.05))
    @test all(get_helix_ends(a, id_a).stop .≈ (0.1, 0.06, 0.125))     # half a turn on: at -y

    # the three axes give the same helix, turned cyclically
    x = mk(); create_helix!(x, 0.05, 0.1, 0.1, 0.04, 0.008, 0.03, 3, DirX(); ε_r = 4.0)
    y = mk(); create_helix!(y, 0.1, 0.05, 0.1, 0.04, 0.008, 0.03, 3, DirY(); ε_r = 4.0)
    z = mk(); create_helix!(z, 0.1, 0.1, 0.05, 0.04, 0.008, 0.03, 3, DirZ(); ε_r = 4.0)
    @test x.material[:, :, :, 2] == permutedims(z.material[:, :, :, 2], (3, 1, 2))
    @test y.material[:, :, :, 2] == permutedims(z.material[:, :, :, 2], (2, 3, 1))

    # invalid requests: nothing recorded
    e = mk()
    @test_throws "helix radius R must be positive" create_helix!(e, 0.1, 0.1, 0.05, 0.0, 0.008, 0.03, 3, DirZ())
    @test_throws "wire radius r must be positive"  create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.0, 0.03, 3, DirZ())
    @test_throws "pitch must be positive"          create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.008, -0.03, 3, DirZ())
    @test_throws "turns must be positive"          create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.008, 0.03, 0, DirZ())
    @test_throws "smaller than the helix radius"   create_helix!(e, 0.1, 0.1, 0.05, 0.01, 0.01, 0.03, 3, DirZ())
    @test_throws "turns would touch"               create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.015, 0.03, 3, DirZ())
    @test_throws "outside the domain"              create_helix!(e, 0.1, 0.3, 0.05, 0.04, 0.008, 0.03, 3, DirZ())
    @test isempty(e._history.steps)
    @test_throws "not a helix" (create_sphere!(e, 0.1, 0.1, 0.1, 0.02); get_helix_ends(e, 1))

    # warnings: turns closer than a cell, clipped, too thin
    @test_logs (:warn, r"can merge") create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.0145, 0.03, 3, DirZ())
    @test_logs (:warn, r"clipped") create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.008, 0.03, 5, DirZ())
    @test_logs (:warn, r"no cells were filled") create_helix!(e, 0.1, 0.1, 0.05, 0.04, 0.0004, 0.03, 3, DirZ())
end
