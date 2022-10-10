@testset "Transforms" begin
  @testset "StdCoords" begin
    trans = StdCoords()
    @test TB.isrevertible(trans)

    # basic tests with Cartesian grid
    grid  = CartesianGrid{T}(10, 10)
    mesh, cache = TB.apply(trans, grid)
    @test all(sides(boundingbox(mesh)) .≤ T(1))
    rgrid = TB.revert(trans, mesh, cache)
    @test rgrid == grid
    mesh2 = TB.reapply(trans, grid, cache)
    @test mesh == mesh2

    # basic tests with views
    vset = view(PointSet(rand(P2, 100)), 1:50)
    vnew, cache = TB.apply(trans, vset)
    @test all(sides(boundingbox(vnew)) .≤ T(1))
    vini = TB.revert(trans, vnew, cache)
    @test vini == vset
    vnew2 = TB.reapply(trans, vset, cache)
    @test vnew == vnew2
  end

  @testset "TaubinSmoothing" begin
    trans = TaubinSmoothing(30)
    @test TB.isrevertible(trans)

    # smoothing doesn't change the topology
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    smesh = trans(mesh)
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)
  end
end
