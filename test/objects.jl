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
