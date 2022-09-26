@testset "Transforms" begin
  @testset "StdCoords" begin
    grid  = CartesianGrid{T}(10, 10)
    trans = StdCoords()
    @test TAPI.isrevertible(trans)

    # basic tests with Cartesian grid
    mesh, cache = TAPI.apply(trans, grid)
    @test all(sides(boundingbox(mesh)) .≤ T(1))
    rgrid = TAPI.revert(trans, mesh, cache)
    @test vertices(rgrid) == vertices(grid)
    mesh2 = TAPI.reapply(trans, grid, cache)
    @test mesh == mesh2
  end

  @testset "TaubinSmoothing" begin
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    trans = TaubinSmoothing(30)
    @test TAPI.isrevertible(trans)

    # smoothing doesn't change the topology
    smesh = trans(mesh)
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)
  end
end
