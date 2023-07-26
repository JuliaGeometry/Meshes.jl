@testset "Transforms" begin
  @testset "Rotate" begin
    @test TB.isrevertible(Rotate)
    @test TB.isinvertible(Rotate)
    @test inv(Rotate(Angle2d(T(π / 2)))) == Rotate(Angle2d(-T(π / 2)))

    # ------
    # POINT
    # ------

    f = Rotate(Angle2d(T(π / 2)))
    g = P2(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ P2(0, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Rotate(Angle2d(T(π / 2)))
    g = Segment(P2(0, 0), P2(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(P2(0, 0), P2(0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Rotate(Angle2d(T(π / 2)))
    g = Box(P2(0, 0), P2(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Quadrangle
    @test r ≈ Quadrangle(P2(0, 0), P2(0, 1), P2(-1, 1), P2(-1, 0))
    q = TB.revert(f, r, c)
    @test q isa Quadrangle
    @test q ≈ convert(Quadrangle, g)

    f = Rotate(V3(1, 0, 0), V3(0, 1, 0))
    g = Box(P3(0, 0, 0), P3(1, 1, 1))
    r, c = TB.apply(f, g)
    @test r isa Hexahedron
    @test r ≈ Hexahedron(
      P3(0, 0, 0),
      P3(0, 1, 0),
      P3(-1, 1, 0),
      P3(-1, 0, 0),
      P3(0, 0, 1),
      P3(0, 1, 1),
      P3(-1, 1, 1),
      P3(-1, 0, 1)
    )
    h = TB.revert(f, r, c)
    @test h isa Hexahedron
    @test h ≈ convert(Hexahedron, g)

    # ----------
    # ROPE/RING
    # ----------

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

    f = Rotate(Angle2d(T(π / 2)))
    p = PolyArea(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    r, c = TB.apply(f, p)
    @test r ≈ PolyArea(P2(0, 0), P2(0, 1), P2(-1, 1), P2(-1, 0))
    @test TB.revert(f, r, c) ≈ p

    # ----------
    # MULTIGEOM
    # ----------

    f = Rotate(Angle2d(T(π / 2)))
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

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

    # ---------
    # POINTSET
    # ---------

    f = Rotate(Angle2d(T(π / 2)))
    d = PointSet([P2(0, 0), P2(1, 0), P2(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([P2(0, 0), P2(0, 1), P2(-1, 1)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Rotate(Angle2d(T(π / 2)))
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d

    # --------------
    # CARTESIANGRID
    # --------------

    f = Rotate(Angle2d(T(π / 2)))
    d = CartesianGrid{T}(10, 10)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Rotate(Angle2d(T(π / 2)))
    p = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d
  end

  @testset "Translate" begin
    @test TB.isrevertible(Translate)
    @test TB.isinvertible(Translate)
    @test inv(Translate(T(1), T(2))) == Translate(T(-1), T(-2))

    # ------
    # POINT
    # ------

    f = Translate(T(1), T(1))
    g = P2(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ P2(2, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Translate(T(1), T(1))
    g = Segment(P2(0, 0), P2(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(P2(1, 1), P2(2, 1))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Translate(T(1), T(1))
    g = Box(P2(0, 0), P2(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(P2(1, 1), P2(2, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Translate(T(1), T(2), T(3))
    g = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(P3(1, 2, 3), P3(2, 2, 3), P3(1, 3, 4))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Translate(T(1), T(1))
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
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

    # ---------
    # POINTSET
    # ---------

    f = Translate(T(1), T(1))
    d = PointSet([P2(0, 0), P2(1, 0), P2(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([P2(1, 1), P2(2, 1), P2(2, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Translate(T(1), T(1))
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d

    # --------------
    # CARTESIANGRID
    # --------------

    f = Translate(T(1), T(1))
    d = CartesianGrid{T}(10, 10)
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(P2(1, 1), P2(11, 11), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Translate(T(1), T(1))
    p = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d
  end

  @testset "Stretch" begin
    @test TB.isrevertible(Stretch)
    @test TB.isinvertible(Stretch)
    @test inv(Stretch(T(1), T(2))) == Stretch(T(1), T(1 / 2))

    # ------
    # POINT
    # ------

    f = Stretch(T(1), T(2))
    g = P2(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ P2(1, 2)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Stretch(T(1), T(2))
    g = Segment(P2(0, 0), P2(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(P2(0, 0), P2(1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Stretch(T(2), T(1))
    g = Segment(P2(0, 0), P2(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(P2(0, 0), P2(2, 0))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Stretch(T(1), T(2))
    g = Box(P2(0, 0), P2(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(P2(0, 0), P2(1, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Stretch(T(1), T(2), T(3))
    g = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 2, 3))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------
    f = Stretch(T(1), T(2))
    f = Translate(T(1), T(1))
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
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

    # ---------
    # POINTSET
    # ---------

    f = Stretch(T(1), T(2))
    d = PointSet([P2(0, 0), P2(1, 0), P2(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([P2(0, 0), P2(1, 0), P2(1, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Stretch(T(1), T(2))
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d

    # --------------
    # CARTESIANGRID
    # --------------

    f = Stretch(T(1), T(2))
    d = CartesianGrid(P2(1, 1), P2(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(P2(1, 2), P2(11, 22), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Stretch(T(1), T(2))
    p = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d
  end

  @testset "StdCoords" begin
    @test TB.isrevertible(StdCoords)

    # ---------
    # POINTSET
    # ---------

    f = StdCoords()
    d = view(PointSet(rand(P2, 100)), 1:50)
    r, c = TB.apply(f, d)
    @test all(sides(boundingbox(r)) .≤ T(1))
    @test TB.revert(f, r, c) ≈ d
    r2 = TB.reapply(f, d, c)
    @test r == r2

    # --------------
    # CARTESIANGRID
    # --------------

    f = StdCoords()
    d = CartesianGrid(P2(1, 1), P2(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(P2(-0.5, -0.5), P2(0.5, 0.5), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    f = StdCoords()
    d = CartesianGrid{T}(10, 20)
    r, c = TB.apply(f, d)
    @test r ≈ CartesianGrid(P2(-0.5, -0.5), P2(0.5, 0.5), dims=(10, 20))
    @test TB.revert(f, r, c) ≈ d
    r2 = TB.reapply(f, d, c)
    @test r == r2
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
