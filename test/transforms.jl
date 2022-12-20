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

  @testset "Smoothing" begin
    # smoothing doesn't change the topology
    trans = LaplaceSmoothing(30)
    @test TB.isrevertible(trans)
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    smesh = trans(mesh)
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)

    # smoothing doesn't change the topology
    trans = TaubinSmoothing(30)
    @test TB.isrevertible(trans)
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    smesh = trans(mesh)
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)
  end

  @testset "Rotation" begin
    # a rotation doesn't change the topology
    trans = Rotation(EulerAngleAxis(T(pi/4), T[1, 0, 0]))
    @test TB.isrevertible(trans)
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    rmesh = trans(mesh)
    @test nvertices(rmesh) == nvertices(mesh)
    @test nelements(rmesh) == nelements(mesh)
    @test topology(rmesh) == topology(mesh)

    # check rotation on a triangle
    triangle = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    trans = Rotation(EulerAngleAxis(T(pi/2), T[0, 0, 1]))
    rtriangle = trans(triangle)
    rpoints = vertices(rtriangle)
    @test rpoints[1] == P3(0, 0, 0)
    @test rpoints[2] == P3(0, 1, 0)
    @test rpoints[3] == P3(-1, 0, 0)
  end

end
