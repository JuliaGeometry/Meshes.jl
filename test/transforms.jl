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
    trans = Rotation(EulerAngleAxis(pi/4, [1, 0, 0]))
    @test TB.isrevertible(trans)
    mesh  = readply(T, joinpath(datadir,"beethoven.ply"))
    rmesh = trans(mesh)
    @test nvertices(rmesh) == nvertices(mesh)
    @test nelements(rmesh) == nelements(mesh)
    @test topology(rmesh) == topology(mesh)

    # check rotation on a triangle
    triangle = Triangle(Point(0,0,0), Point(1,0,0), Point(0,1,0))
    trans = Rotation(EulerAngleAxis(-pi/2, [0, 0, 1]))
    rtriangle = trans(triangle)
    @test in(Point(0,0,0), rtriangle)
    @test in(Point(-0.01,0.9,0), rtriangle)
    @test in(Point(-0.9,0.01,0), rtriangle)
  end

end
