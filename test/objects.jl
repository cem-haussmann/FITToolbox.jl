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

@testset "create_brick!: reaching beyond the domain" begin
    # clipped there with a warning, like the round solids, and recorded as given
    d = create_domain([1.0, 1.0, 1.0], [0.1, 0.1, 0.1])
    id = @test_logs (:warn, r"clipped") create_brick!(d, 0.6, 0.0, 0.0, 0.8, 1.0, 0.2;
                                                      ε_r = 2.0, name = "plate")
    @test id == 1 && only(list_objects(d)).object.lengths == (0.8, 1.0, 0.2)
    filled = d.material[:, :, :, 2] .== 2.0
    @test all(filled[7:10, :, 1:2]) && count(filled) == 4 * 10 * 2   # u = 0.6 … 1.0 m, w = 0 … 0.2 m
    # clipped again on another grid, warning once
    fine = @test_logs (:warn, r"clipped") change_resolution(d, 0.05)
    @test count(fine.material[:, :, :, 2] .== 2.0) == 8 * 20 * 4

    # still an error, nothing recorded: the corner outside, a length that is not
    # positive, or less than half a cell left inside the domain (no warning then)
    @test_throws "corner (1.2, 0.0, 0.0) m lies outside" create_brick!(d, 1.2, 0.0, 0.0, 0.2, 0.2, 0.2)
    @test_throws "lies outside the domain" create_brick!(d, 0.5, -0.1, 0.0, 0.2, 0.5, 0.2)
    @test_throws "lengths must be positive" create_brick!(d, 0.5, 0.5, 0.5, -0.2, 0.2, 0.2)
    @test_throws "lengths must be positive" create_brick!(d, 0.5, 0.5, 0.5, 0.2, 0.0, 0.2)
    @test_logs @test_throws "zero extent in u" create_brick!(d, 0.98, 0.0, 0.0, 0.5, 0.2, 0.2)
    @test_throws "zero extent in w" create_brick!(d, 0.5, 0.5, 0.5, 0.2, 0.2, 0.03)
    @test length(d._history.steps) == 1
end

@testset "create_brick!: snaps like the loops and charges" begin
    # a face exactly halfway between two nodes goes to the lower one, as the corners of a
    # rectangular loop and a charge do: the same coordinates give the same nodes (quarter
    # metre cells, so that the halfway points are exact)
    d = create_domain([2.0, 2.0, 2.0], [0.25, 0.25, 0.25])
    create_brick!(d, 0.125, 0.375, 0.0, 0.5, 0.25, 2.0; ε_r = 2.0)
    filled = d.material[:, :, 1, 2] .== 2.0
    @test findall(vec(any(filled; dims = 2))) == [1, 2]                # u: 0 … 0.5, not 0.25 … 0.75
    @test findall(vec(any(filled; dims = 1))) == [2]                   # v: 0.25 … 0.5
    J = get_source(d, create_rectangular_loop_source!(d, 0.125, 0.375, 1.0, 0.5, 0.25, DirZ()))
    starts = [get_position_of_index(d, PrimalNode(), mod1(p, d.Np)) for p in findall(!iszero, J)]
    @test extrema(s[1] for s in starts) == (0.0, 0.5) && extrema(s[2] for s in starts) == (0.25, 0.5)
    q = get_source(d, create_electric_charge!(d, 0.125, 0.375, 1.0; charge = 1.0))
    @test get_position_of_index(d, PrimalNode(), findfirst(!iszero, q))[1:2] == (0.0, 0.25)
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

@testset "the round solids scan only the cells around them" begin
    # the cells that can hold part of a solid within [lo, hi]: those overlapping it and
    # one more on each side, within the grid
    nodes = [0.0, 0.1, 0.3, 0.4, 0.7, 1.0]                             # five graded cells
    @test FITToolbox._cell_range(nodes, 0.35, 0.5) == 2:5             # cells 3 and 4, ± 1
    @test FITToolbox._cell_range(nodes, 0.3, 0.4) == 2:5              # on nodes: cell 3, both sides
    @test FITToolbox._cell_range(nodes, -2.0, 0.05) == 1:3            # reaching below the grid
    @test FITToolbox._cell_range(nodes, 0.95, 3.0) == 4:5             # and above it
    @test FITToolbox._cell_range(nodes, 2.0, 3.0) == 5:5              # beyond: nothing to fill

    # The same material as scanning the whole grid: every cell outside the range has none
    # of its 27 samples inside. A sphere and a helix, in the middle and clipped at a face.
    d = create_domain([0.1, 0.05, 0.2, 0.05, 0.1, 0.15, 0.1], [0.08, 0.12, 0.1, 0.05, 0.15],
                      [0.1, 0.1, 0.05, 0.15, 0.1, 0.1])
    solids = [FITToolbox.Sphere((0.37, 0.21, 0.3), 0.17, FITToolbox.Material(1.0, 1.0, 1.0)),
              FITToolbox.Sphere((0.05, 0.4, 0.55), 0.2, FITToolbox.Material(1.0, 1.0, 1.0)),
              FITToolbox.Helix((0.4, 0.25, 0.1), 0.12, 0.03, 0.1, 3.5, DirZ(), 30.0, true,
                               FITToolbox.Material(1.0, 1.0, 1.0))]
    for s in solids
        lo, hi = FITToolbox._extent(s)
        inside = FITToolbox._inside_test(s)
        ref = FITToolbox._reference(s)
        cells = FITToolbox._cell_ranges(d, (lo, hi))
        nodes = (d.nodes_u, d.nodes_v, d.nodes_w)
        hit(c) = any(Iterators.product(ntuple(a -> [nodes[a][c[a]] + f * (nodes[a][c[a]+1] - nodes[a][c[a]]) - ref[a]
                                                     for f in (1/6, 1/2, 5/6)], 3)...)) do p
            inside(p)
        end
        outside = [c for c in CartesianIndices((d.Nu - 1, d.Nv - 1, d.Nw - 1))
                   if !all(a -> c[a] in cells[a], 1:3)]
        @test !isempty(outside) && !any(c -> hit(Tuple(c)), outside)
    end
end

@testset "objects reaching the boundary exactly" begin
    # ten cells of 0.1 m end at 0.9999999999999999 m: objects up to 1.0 m must still fit
    d = create_domain([1.0, 1.0, 1.0], [0.1, 0.1, 0.1])
    @test d.nodes_u[end] < 1.0                                         # the rounding this is about
    @test (@test_logs create_brick!(d, 0.0, 0.0, 0.0, 1.0, 1.0, 0.1; ε_r = 2.0)) == 1   # not clipped
    @test all(d.material[:, :, 1, 2] .== 2.0)                         # every cell of the bottom layer
    @test_logs create_cylinder!(d, 0.5, 0.5, 0.0, 0.2, 1.0, DirZ())   # up to the top: not clipped
    @test_logs (:warn, r"clipped") create_sphere!(d, 1.0, 0.5, 0.5, 0.2)   # centre on the boundary
    @test create_electric_charge!(d, 1.0, 1.0, 1.0; charge = 1e-9) == 4
    @test get_source(d, 4)[end] == 1e-9                               # the last node

    # beyond it: a brick is clipped, with a warning; a charge is an error
    @test_logs (:warn, r"clipped") create_brick!(d, 0.0, 0.0, 0.0, 1.001, 1.0, 0.1)
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

@testset "signed distances of the solids" begin
    # agree with the inside test in sign, are distances (|∇d| = 1 near the surface), and
    # take known values at the centre lines (cells smaller than the thinnest wire)
    d = create_domain([1.0, 1.0, 1.0], [0.025, 0.025, 0.025])
    create_brick!(d, 0.2, 0.3, 0.1, 0.4, 0.2, 0.5)
    create_sphere!(d, 0.5, 0.5, 0.5, 0.2)
    create_cylinder!(d, 0.5, 0.5, 0.2, 0.1, 0.4, DirY())
    create_torus!(d, 0.5, 0.5, 0.5, 0.25, 0.06, DirX())
    create_helix!(d, 0.5, 0.5, 0.2, 0.2, 0.03, 0.2, 3.3, DirZ(); start_angle = 40, left_handed = true)
    grad(o, p; h = 1e-6) = ntuple(k -> (FITToolbox._signed_distance(o, p .+ ntuple(i -> i == k ? h : 0.0, 3)) -
                                        FITToolbox._signed_distance(o, p)) / h, 3)
    for e in list_objects(d).entries
        o = e.object
        ref, inside = FITToolbox._reference(o), FITToolbox._inside_test(o)
        lo, hi = FITToolbox._extent(o)
        disagreements, grads = 0, Float64[]
        # a lattice of 17³ points over the extent and a margin, off the grid planes
        for z in range(lo[3] - 0.05, hi[3] + 0.05; length = 17),
            y in range(lo[2] - 0.05, hi[2] + 0.05; length = 17),
            x in range(lo[1] - 0.05, hi[1] + 0.05; length = 17)
            p = (x, y, z) .- ref .+ 1e-3
            sd = FITToolbox._signed_distance(o, p)
            (sd <= 0) == inside(p) || (disagreements += 1)
            abs(sd) < 0.02 && push!(grads, sqrt(sum(abs2, grad(o, p))))
        end
        @test disagreements == 0
        @test !isempty(grads) && all(g -> isapprox(g, 1.0; atol = 1e-3), grads)
        # outside the extent is outside the solid
        @test FITToolbox._signed_distance(o, hi .+ 0.1 .- ref) > 0
    end
    brick, sphere, cyl, torus, helix = (e.object for e in list_objects(d).entries)
    @test FITToolbox._signed_distance(brick, (0.2, 0.1, 0.25)) ≈ -0.1          # the centre: half the shortest side
    @test FITToolbox._signed_distance(sphere, (0.0, 0.0, 0.0)) ≈ -0.2
    @test FITToolbox._signed_distance(cyl, (0.0, 0.2, 0.0)) ≈ -0.1             # on the axis, mid-height: the radius
    @test FITToolbox._signed_distance(cyl, (0.0, 0.05, 0.0)) ≈ -0.05           # near the base: the base is closer
    @test FITToolbox._signed_distance(torus, (0.0, 0.25, 0.0)) ≈ -0.06         # on the tube's centre circle
    @test FITToolbox._signed_distance(torus, (0.0, 0.0, 0.0)) ≈ 0.25 - 0.06    # the ring's centre: in the hole
    on_line = FITToolbox._helix_point(helix, 2.0) .- helix.base
    @test FITToolbox._signed_distance(helix, on_line) ≈ -0.03
    @test FITToolbox._signed_distance(helix, (0.0, 0.0, 0.1)) ≈ 0.2 - 0.03    # on the axis: R - r away
end
