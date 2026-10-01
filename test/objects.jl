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
