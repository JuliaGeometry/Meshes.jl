@testset "Transforms" begin
  @testset "Rotate" begin
    # check rotation on a triangle
    tri  = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    rtri = tri |> Rotate(EulerAngleAxis(T(pi/2), T[0, 0, 1]))
    rpts = vertices(rtri)
    @test rpts[1] ≈ P3(0, 0, 0)
    @test rpts[2] ≈ P3(0, 1, 0)
    @test rpts[3] ≈ P3(-1, 0, 0)
    
    # triangle in the plane z=0
    tri  = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    # rotate to the plane x=0
    rtri = tri |> Rotate(V3(0, 0, 1), V3(1, 0, 0))
    # check that the rotated triangle is in the x=0 plane
    rpts = coordinates.(vertices(rtri))
    @test isapprox(rpts[1][1], zero(T); atol = atol(T))
    @test isapprox(rpts[2][1], zero(T); atol = atol(T))
    @test isapprox(rpts[3][1], zero(T); atol = atol(T))
  end

  @testset "Translate" begin
    # check translation on a triangle
    tri  = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    ttri = tri |> Translate(T(1), T(2), T(3))
    tpts = vertices(ttri)
    @test tpts[1] ≈ P3(1, 2, 3)
    @test tpts[2] ≈ P3(2, 2, 3)
    @test tpts[3] ≈ P3(1, 3, 4)
  end

  @testset "Stretch" begin
    # check scaling on a triangle
    tri  = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    stri = tri |> Stretch(T(1), T(2), T(3))
    spts = vertices(stri)
    @test spts[1] ≈ P3(0, 0, 0)
    @test spts[2] ≈ P3(1, 0, 0)
    @test spts[3] ≈ P3(0, 2, 3)
  end
  
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
    @test all(vini .≈ vset)
    vnew2 = TB.reapply(trans, vset, cache)
    @test vnew == vnew2
  end

  @testset "Repair{1}" begin
    # a tetrahedron with an unused vertex
    verts = P3[(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)]
    connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)], Triangle)
    mesh  = SimpleMesh(verts, connec)
    rmesh = mesh |> Repair{1}()
    @test nvertices(rmesh) == nvertices(mesh) - 1
    @test nelements(rmesh) == nelements(mesh)
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
end
