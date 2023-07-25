@testset "Transforms" begin
  @testset "Rotate" begin
    # ------
    # CHAIN
    # ------

    f = Rotate(Angle2d(T(π / 2)))
    g = Rope(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Rope(P2(0, 0), P2(0, 1), P2(-1, 1), P2(-1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Rotate(Angle2d(T(π / 2)))
    g = Ring(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Ring(P2(0, 0), P2(0, 1), P2(-1, 1), P2(-1, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Rotate(AngleAxis(T(π / 2), T(0), T(0), T(1)))
    g = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(P3(0, 0, 0), P3(0, 1, 0), P3(-1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Rotate(V3(0, 0, 1), V3(1, 0, 0))
    g = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(P3(0, 0, 0), P3(0, 0, -1), P3(0, 1, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POLYAREA
    # ---------

    #f = Rotate(Angle2d(T(π / 2)))
    #p = PolyArea(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    #r, c = TB.apply(f, p)
    #@test r ≈ PolyArea(P2(0, 0), P2(0, 1), P2(-1, 1), P2(-1, 0))
    #@test TB.revert(f, r, c) ≈ p

    # ------
    # PLANE
    # ------

    f = Rotate(V3(0, 0, 1), V3(1, 0, 0))
    g = Plane(P3(0, 0, 0), V3(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(P3(0, 0, 0), V3(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Rotate(V3(0, 0, 1), V3(1, 0, 0))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(P3(0, 0, 0), P3(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g
  end

  @testset "Translate" begin
    # ---------
    # TRIANGLE
    # ---------

    f = Translate(T(1), T(2), T(3))
    g = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(P3(1, 2, 3), P3(2, 2, 3), P3(1, 3, 4))
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Translate(T(0), T(0), T(1))
    g = Plane(P3(0, 0, 0), V3(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(P3(0, 0, 0), V3(0, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Translate(T(0), T(0), T(1))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(P3(0, 0, 1), P3(0, 0, 2))
    @test TB.revert(f, r, c) ≈ g
  end

  @testset "Stretch" begin
    # ---------
    # TRIANGLE
    # ---------

    f = Stretch(T(1), T(2), T(3))
    g = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 2, 3))
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Stretch(T(1), T(1), T(2))
    g = Plane(P3(1, 1, 1), V3(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    f = Stretch(T(2), T(1), T(1))
    g = Plane(P3(1, 1, 1), V3(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Stretch(T(1), T(1), T(2))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(P3(0, 0, 1), P3(0, 0, 2))
    @test TB.revert(f, r, c) ≈ g
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
