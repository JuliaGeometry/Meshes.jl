@testset "Transforms" begin
  @testset "Rotate" begin
    # rotate around z axis
    f = Rotate(AngleAxis(T(pi / 2), T(0), T(0), T(1)))
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    r, c = TB.apply(f, t)
    @test r ≈ Triangle(P3(0, 0, 0), P3(0, 1, 0), P3(-1, 0, 0))
    @test TB.revert(f, r, c) ≈ t

    # rotate from plane z=0 to plane x=0
    f = Rotate(V3(0, 0, 1), V3(1, 0, 0))
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    r, c = TB.apply(f, t)
    @test r ≈ Triangle(P3(0, 0, 0), P3(0, 0, -1), P3(0, 1, 0))
    @test TB.revert(f, r, c) ≈ t

    # rotate plane
    f = Rotate(V3(0, 0, 1), V3(1, 0, 0))
    p = Plane(P3(0, 0, 0), V3(0, 0, 1))
    r, c = TB.apply(f, p)
    @test r ≈ Plane(P3(0, 0, 0), V3(1, 0, 0))
    @test TB.revert(f, r, c) ≈ p
  end

  @testset "Translate" begin
    tri = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    ttri = tri |> Translate(T(1), T(2), T(3))
    tpts = vertices(ttri)
    @test tpts[1] ≈ P3(1, 2, 3)
    @test tpts[2] ≈ P3(2, 2, 3)
    @test tpts[3] ≈ P3(1, 3, 4)
  end

  @testset "Stretch" begin
    # check scaling on a triangle
    tri = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
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
    grid = CartesianGrid{T}(10, 10)
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

  @testset "Repair{0}" begin
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 0), (1, 1), (0, 1), (0, 1)])
    rpoly = poly |> Repair{0}()
    @test nvertices(rpoly) == 4
    @test vertices(rpoly) == P2[(0, 0), (1, 0), (1, 1), (0, 1)]
  end

  @testset "Repair{1}" begin
    # a tetrahedron with an unused vertex
    points = P3[(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)]
    connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)])
    mesh = SimpleMesh(points, connec)
    rmesh = mesh |> Repair{1}()
    @test nvertices(rmesh) == nvertices(mesh) - 1
    @test nelements(rmesh) == nelements(mesh)
    @test P3(5, 5, 5) ∉ vertices(rmesh)
  end

  @testset "Repair{7}" begin
    # mesh with inconsistent orientation
    points = rand(P3, 6)
    connec = connect.([(1, 2, 3), (3, 4, 2), (4, 3, 5), (6, 3, 1)])
    mesh = SimpleMesh(points, connec)
    rmesh = mesh |> Repair{7}()
    topo = topology(mesh)
    rtopo = topology(rmesh)
    e = collect(elements(topo))
    n = collect(elements(rtopo))
    @test n[1] == e[1]
    @test n[2] != e[2]
    @test n[3] != e[3]
    @test n[4] != e[4]
  end

  @testset "Repair{8}" begin
    poly =
      PolyArea(P2[(0.0, 0.0), (0.5, -0.5), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (-0.5, 0.5)])
    rpoly = poly |> Repair{8}()
    @test nvertices(rpoly) == 4
    @test vertices(rpoly) == P2[(0.5, -0.5), (1.5, 0.5), (0.5, 1.5), (-0.5, 0.5)]
  end

  @testset "Smoothing" begin
    # smoothing doesn't change the topology
    trans = LaplaceSmoothing(30)
    @test TB.isrevertible(trans)
    mesh = readply(T, joinpath(datadir, "beethoven.ply"))
    smesh = trans(mesh)
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)

    # smoothing doesn't change the topology
    trans = TaubinSmoothing(30)
    @test TB.isrevertible(trans)
    mesh = readply(T, joinpath(datadir, "beethoven.ply"))
    smesh = trans(mesh)
    @test nvertices(smesh) == nvertices(mesh)
    @test nelements(smesh) == nelements(mesh)
    @test topology(smesh) == topology(mesh)
  end
end
