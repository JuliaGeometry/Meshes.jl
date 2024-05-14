@testset "Transforms" begin
  @testset "Rotate" begin
    @test isaffine(Rotate)
    @test TB.isrevertible(Rotate)
    @test TB.isinvertible(Rotate)
    @test TB.inverse(Rotate(Angle2d(T(π / 2)))) == Rotate(Angle2d(-T(π / 2)))
    rot = Angle2d(T(π / 2))
    f = Rotate(rot)
    @test TB.parameters(f) == (; rot)

    # ----
    # VEC
    # ----

    f = Rotate(Angle2d(T(π / 2)))
    v = vector(1, 0)
    r, c = TB.apply(f, v)
    @test r ≈ vector(0, 1)
    @test TB.revert(f, r, c) ≈ v

    # ------
    # POINT
    # ------

    f = Rotate(Angle2d(T(π / 2)))
    g = point(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ point(0, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Rotate(Angle2d(T(π / 2)))
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(0, 0), point(0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Rotate(Angle2d(T(π / 2)))
    g = Box(point(0, 0), point(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Quadrangle
    @test r ≈ Quadrangle(point(0, 0), point(0, 1), point(-1, 1), point(-1, 0))
    q = TB.revert(f, r, c)
    @test q isa Quadrangle
    @test q ≈ convert(Quadrangle, g)

    f = Rotate(vector(1, 0, 0), vector(0, 1, 0))
    g = Box(point(0, 0, 0), point(1, 1, 1))
    r, c = TB.apply(f, g)
    @test r isa Hexahedron
    @test r ≈ Hexahedron(
      point(0, 0, 0),
      point(0, 1, 0),
      point(-1, 1, 0),
      point(-1, 0, 0),
      point(0, 0, 1),
      point(0, 1, 1),
      point(-1, 1, 1),
      point(-1, 0, 1)
    )
    h = TB.revert(f, r, c)
    @test h isa Hexahedron
    @test h ≈ convert(Hexahedron, g)

    # ----------
    # ROPE/RING
    # ----------

    f = Rotate(Angle2d(T(π / 2)))
    g = Rope(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Rope(point(0, 0), point(0, 1), point(-1, 1), point(-1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Rotate(Angle2d(T(π / 2)))
    g = Ring(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Ring(point(0, 0), point(0, 1), point(-1, 1), point(-1, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Rotate(AngleAxis(T(π / 2), T(0), T(0), T(1)))
    g = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(point(0, 0, 0), point(0, 1, 0), point(-1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Rotate(vector(0, 0, 1), vector(1, 0, 0))
    g = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(point(0, 0, 0), point(0, 0, -1), point(0, 1, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POLYAREA
    # ---------

    f = Rotate(Angle2d(T(π / 2)))
    p = PolyArea(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    r, c = TB.apply(f, p)
    @test r ≈ PolyArea(point(0, 0), point(0, 1), point(-1, 1), point(-1, 0))
    @test TB.revert(f, r, c) ≈ p

    # ----------
    # MULTIGEOM
    # ----------

    f = Rotate(Angle2d(T(π / 2)))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Rotate(vector(0, 0, 1), vector(1, 0, 0))
    g = Plane(point(0, 0, 0), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(point(0, 0, 0), vector(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Rotate(vector(0, 0, 1), vector(1, 0, 0))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(point(0, 0, 0), point(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Rotate(Angle2d(T(π / 2)))
    d = PointSet([point(0, 0), point(1, 0), point(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([point(0, 0), point(0, 1), point(-1, 1)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Rotate(Angle2d(T(π / 2)))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])
    @test all(TB.revert(f, r, c) .≈ d)

    # --------------
    # CARTESIANGRID
    # --------------

    f = Rotate(Angle2d(T(π / 2)))
    d = cartgrid(10, 10)
    r, c = TB.apply(f, d)
    @test r isa TransformedGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Rotate(Angle2d(T(π / 2)))
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa TransformedGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Rotate(Angle2d(T(π / 2)))
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa TransformedGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Rotate(Angle2d(T(π / 2)))
    p = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ---------
    # FALLBACK
    # ---------

    f = Rotate(T(π / 2))
    v = vector(1, 0)
    r, c = TB.apply(f, v)
    @test r ≈ vector(0, 1)
    @test TB.revert(f, r, c) ≈ v
  end

  @testset "Translate" begin
    @test isaffine(Translate)
    @test TB.isrevertible(Translate)
    @test TB.isinvertible(Translate)
    @test TB.inverse(Translate(T(1), T(2))) == Translate(T(-1), T(-2))
    offsets = (T(1) * u"m", T(2) * u"m")
    f = Translate(offsets)
    @test TB.parameters(f) == (; offsets)
    f = Translate(T(1), T(2))
    @test TB.parameters(f) == (; offsets)
    f = Translate(T(1), 2)
    @test TB.parameters(f) == (; offsets)
    f = Translate(1, 2)
    @test TB.parameters(f) == (; offsets)

    # ----
    # VEC
    # ----

    f = Translate(T(1), T(1))
    v = vector(1, 0)
    r, c = TB.apply(f, v)
    @test r ≈ vector(1, 0)
    @test TB.revert(f, r, c) ≈ v

    # ------
    # POINT
    # ------

    f = Translate(T(1), T(1))
    g = point(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ point(2, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Translate(T(1), T(1))
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(1, 1), point(2, 1))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Translate(T(1), T(1))
    g = Box(point(0, 0), point(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(point(1, 1), point(2, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Translate(T(1), T(2), T(3))
    g = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(point(1, 2, 3), point(2, 2, 3), point(1, 3, 4))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Translate(T(1), T(1))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Translate(T(0), T(0), T(1))
    g = Plane(point(0, 0, 0), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(point(0, 0, 1), vector(0, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Translate(T(0), T(0), T(1))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(point(0, 0, 1), point(0, 0, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Translate(T(1), T(1))
    d = PointSet([point(0, 0), point(1, 0), point(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([point(1, 1), point(2, 1), point(2, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Translate(T(1), T(1))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])
    @test all(TB.revert(f, r, c) .≈ d)

    # --------------
    # CARTESIANGRID
    # --------------

    f = Translate(T(1), T(1))
    d = cartgrid(10, 10)
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(point(1, 1), point(11, 11), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Translate(T(1), T(1))
    d = RectilinearGrid(T.(0:10), T.(0:10))
    r, c = TB.apply(f, d)
    @test r isa RectilinearGrid
    @test r ≈ RectilinearGrid(T.(1:11), T.(1:11))
    @test TB.revert(f, r, c) ≈ d

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Translate(T(1), T(1))
    d = StructuredGrid(repeat(T.(0:10), 1, 11), repeat(T.(0:10)', 11, 1))
    r, c = TB.apply(f, d)
    @test r isa StructuredGrid
    @test r ≈ StructuredGrid(repeat(T.(1:11), 1, 11), repeat(T.(1:11)', 11, 1))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Translate(T(1), T(1))
    p = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d
  end

  @testset "Affine" begin
    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    @test isaffine(f)
    @test TB.isrevertible(f)
    @test TB.isinvertible(f)
    @test TB.inverse(f) == Affine(Angle2d(-T(π / 2)), Angle2d(-T(π / 2)) * -T[1, 1])
    f = Affine(T[6 3; 10 5], T[1, 1])
    @test !TB.isrevertible(f)
    @test !TB.isinvertible(f)
    A, b = Angle2d(T(π / 2)), SVector(T(1) * u"m", T(1) * u"m")
    f = Affine(A, b)
    @test TB.parameters(f) == (; A, b)
    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    @test TB.parameters(f) == (; A, b)
    f = Affine(Angle2d(T(π / 2)), [1, 1])
    @test TB.parameters(f) == (; A, b)

    # ----
    # VEC
    # ----

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    v = vector(1, 0)
    r, c = TB.apply(f, v)
    @test r ≈ vector(0, 1)
    @test TB.revert(f, r, c) ≈ v

    # ------
    # POINT
    # ------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    g = point(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ point(1, 2)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(1, 1), point(1, 2))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    g = Box(point(0, 0), point(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Quadrangle
    @test r ≈ Quadrangle(point(1, 1), point(1, 2), point(0, 2), point(0, 1))
    q = TB.revert(f, r, c)
    @test q isa Quadrangle
    @test q ≈ convert(Quadrangle, g)

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[1, 2, 3])
    g = Box(point(0, 0, 0), point(1, 1, 1))
    r, c = TB.apply(f, g)
    @test r isa Hexahedron
    @test r ≈ Hexahedron(
      point(1, 2, 3),
      point(1, 2, 2),
      point(1, 3, 2),
      point(1, 3, 3),
      point(2, 2, 3),
      point(2, 2, 2),
      point(2, 3, 2),
      point(2, 3, 3)
    )
    h = TB.revert(f, r, c)
    @test h isa Hexahedron
    @test h ≈ convert(Hexahedron, g)

    # ---------
    # TRIANGLE
    # ---------

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[1, 2, 3])
    g = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(point(1, 2, 3), point(1, 2, 2), point(2, 3, 3))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[0, 0, 1])
    g = Plane(point(0, 0, 0), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(point(0, 0, 1), vector(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[0, 0, 1])
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(point(0, 0, 1), point(1, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    d = PointSet([point(0, 0), point(1, 0), point(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([point(1, 1), point(1, 2), point(0, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])
    @test all(TB.revert(f, r, c) .≈ d)

    # --------------
    # CARTESIANGRID
    # --------------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    d = cartgrid(10, 10)
    r, c = TB.apply(f, d)
    @test r isa TransformedGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa TransformedGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa TransformedGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ----------
    # TRANSFORM
    # ----------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    s = Rotate(T(π / 2)) → Translate(T(1), T(1))
    v = vector(1, 0)
    g1 = point(1, 0)
    g2 = Segment(point(0, 0), point(1, 0))
    g3 = Box(point(0, 0), point(1, 1))
    @test f(v) ≈ s(v)
    @test f(g1) ≈ s(g1)
    @test f(g2) ≈ s(g2)
    @test f(g3) ≈ s(g3)

    # ------------
    # CONSTRUCTOR
    # ------------

    # conversion to SArray
    f = Affine(T[0 -1; 1 0], SVector{2}(T[1, 1]))
    @test f.A isa SMatrix
    f = Affine(SMatrix{2,2}(T[0 -1; 1 0]), T[1, 1])
    @test f.b isa SVector
    f = Affine(T[0 -1; 1 0], T[1, 1])
    @test f.A isa SMatrix
    @test f.b isa SVector

    # error: A must be a square matrix
    @test_throws ArgumentError Affine(T[1 1; 2 2; 3 3], T[1, 2])
    # error: A and b must have the same dimension
    @test_throws ArgumentError Affine(T[1 1; 2 2], T[1, 2, 3])
  end

  @testset "Scale" begin
    @test isaffine(Scale)
    @test TB.isrevertible(Scale)
    @test TB.isinvertible(Scale)
    @test TB.inverse(Scale(T(1), T(2))) == Scale(T(1), T(1 / 2))
    factors = (T(1), T(2))
    f = Scale(factors)
    @test TB.parameters(f) == (; factors)

    # ----
    # VEC
    # ----

    f = Scale(T(1), T(2))
    v = vector(1, 1)
    r, c = TB.apply(f, v)
    @test r ≈ vector(1, 2)
    @test TB.revert(f, r, c) ≈ v

    # ------
    # POINT
    # ------

    f = Scale(T(1), T(2))
    g = point(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ point(1, 2)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Scale(T(1), T(2))
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(0, 0), point(1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Scale(T(2), T(1))
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(0, 0), point(2, 0))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Scale(T(1), T(2))
    g = Box(point(0, 0), point(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(point(0, 0), point(1, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Scale(T(1), T(2), T(3))
    g = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 2, 3))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Scale(T(1), T(2))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Scale(T(1), T(1), T(2))
    g = Plane(point(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(point(1, 1, 2), vector(0, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    f = Scale(T(2), T(1), T(1))
    g = Plane(point(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Scale(T(1), T(1), T(2))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(point(0, 0, 0), point(0, 0, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Scale(T(1), T(2))
    d = PointSet([point(0, 0), point(1, 0), point(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([point(0, 0), point(1, 0), point(1, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Scale(T(1), T(2))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])
    @test all(TB.revert(f, r, c) .≈ d)

    # --------------
    # CARTESIANGRID
    # --------------

    f = Scale(T(1), T(2))
    d = CartesianGrid(point(1, 1), point(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(point(1, 2), point(11, 22), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Scale(T(1), T(2))
    p = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d
  end

  @testset "Stretch" begin
    @test !isaffine(Stretch)
    @test TB.isrevertible(Stretch)
    @test TB.isinvertible(Stretch)
    @test TB.inverse(Stretch(T(1), T(2))) == Stretch(T(1), T(1 / 2))
    factors = (T(1), T(2))
    f = Stretch(factors)
    @test TB.parameters(f) == (; factors)

    # ----
    # VEC
    # ----

    f = Stretch(T(1), T(2))
    v = vector(1, 1)
    r, c = TB.apply(f, v)
    @test r ≈ vector(1, 2)
    @test TB.revert(f, r, c) ≈ v

    # ------
    # POINT
    # ------

    f = Stretch(T(1), T(2))
    g = point(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ point(1, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Stretch(T(1), T(2))
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(0, 0), point(1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Stretch(T(2), T(1))
    g = Segment(point(0, 0), point(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(point(-0.5, 0), point(1.5, 0))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Stretch(T(1), T(2))
    g = Box(point(0, 0), point(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(point(0, -0.5), point(1, 1.5))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Stretch(T(1), T(2), T(2))
    g = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(point(0, -1 / 3, -1 / 3), point(1, -1 / 3, -1 / 3), point(0, 10 / 6, 10 / 6))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Stretch(T(1), T(2))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Stretch(T(1), T(1), T(2))
    g = Plane(point(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    f = Stretch(T(2), T(1), T(1))
    g = Plane(point(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Stretch(T(1), T(1), T(2))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(point(0, 0, -0.5), point(0, 0, 1.5))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Stretch(T(1), T(2))
    d = PointSet([point(0, 0), point(1, 0), point(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([point(0, -1 / 3), point(1, -1 / 3), point(1, 10 / 6)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Stretch(T(1), T(2))
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ d
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])
    @test all(TB.revert(f, r, c) .≈ d)

    # --------------
    # CARTESIANGRID
    # --------------

    f = Stretch(T(1), T(2))
    d = CartesianGrid(point(1, 1), point(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(point(1, -4), point(11, 16), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Stretch(T(1), T(2))
    p = point.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d
  end

  @testset "StdCoords" begin
    @test !isaffine(StdCoords)
    @test TB.isrevertible(StdCoords)

    # ---------
    # POINTSET
    # ---------

    f = StdCoords()
    d = view(PointSet(randpoint2(100)), 1:50)
    r, c = TB.apply(f, d)
    @test all(sides(boundingbox(r)) .≤ oneunit(ℳ))
    @test TB.revert(f, r, c) ≈ d
    r2 = TB.reapply(f, d, c)
    @test r == r2

    # --------------
    # CARTESIANGRID
    # --------------

    f = StdCoords()
    d = CartesianGrid(point(1, 1), point(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(point(-0.5, -0.5), point(0.5, 0.5), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    f = StdCoords()
    d = cartgrid(10, 20)
    r, c = TB.apply(f, d)
    @test r ≈ CartesianGrid(point(-0.5, -0.5), point(0.5, 0.5), dims=(10, 20))
    @test TB.revert(f, r, c) ≈ d
    r2 = TB.reapply(f, d, c)
    @test r == r2
  end

  @testset "Repair{0}" begin
    @test !isaffine(Repair)
    poly = PolyArea(point.([(0, 0), (1, 0), (1, 0), (1, 1), (0, 1), (0, 1)]))
    rpoly = poly |> Repair{0}()
    @test nvertices(rpoly) == 4
    @test vertices(rpoly) == point.([(0, 0), (1, 0), (1, 1), (0, 1)])
  end

  @testset "Repair{1}" begin
    # a tetrahedron with an unused vertex
    points = point.([(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)])
    connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)])
    mesh = SimpleMesh(points, connec)
    rmesh = mesh |> Repair{1}()
    @test nvertices(rmesh) == nvertices(mesh) - 1
    @test nelements(rmesh) == nelements(mesh)
    @test point(5, 5, 5) ∉ vertices(rmesh)
  end

  @testset "Repair{2}" begin end

  @testset "Repair{3}" begin end

  @testset "Repair{4}" begin end

  @testset "Repair{5}" begin end

  @testset "Repair{6}" begin end

  @testset "Repair{7}" begin
    # mesh with inconsistent orientation
    points = randpoint3(6)
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
    poly = PolyArea(
      point.([(0.0, 0.0), (0.5, -0.5), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (-0.5, 0.5)])
    )
    rpoly = poly |> Repair{8}()
    @test nvertices(rpoly) == 4
    @test vertices(rpoly) == point.([(0.5, -0.5), (1.5, 0.5), (0.5, 1.5), (-0.5, 0.5)])

    # degenerate triangle with repeated vertices
    poly = PolyArea(point.([(0, 0), (1, 1), (1, 1)]))
    rpoly = poly |> Repair{8}()
    @test !hasholes(rpoly)
    @test rings(rpoly) == [Ring(point(0, 0))]
    @test vertices(rpoly) == [point(0, 0)]
  end

  @testset "Repair{9}" begin
    poly = Quadrangle(point(0, 1, 0), point(1, 1, 0), point(1, 0, 0), point(0, 0, 0))
    bpoly = poly |> Repair{9}()
    @test bpoly isa Quadrangle
    @test bpoly == poly
  end

  @testset "Repair{10}" begin
    outer = Ring(point.([(0, 0), (0, 3), (2, 3), (2, 2), (3, 2), (3, 0)]))
    inner = Ring(point.([(1, 1), (1, 2), (2, 2), (2, 1)]))
    poly = PolyArea(outer, inner)
    repair = Repair{10}()
    rpoly, cache = TB.apply(repair, poly)
    @test nvertices(rpoly) == nvertices(poly)
    @test length(first(rings(rpoly))) > length(first(rings(poly)))
    opoly = TB.revert(repair, rpoly, cache)
    @test opoly == poly
  end

  @testset "Bridge" begin
    @test !isaffine(Bridge)
    δ = T(0.01) * u"m"
    f = Bridge(δ)
    @test TB.parameters(f) == (; δ)
    f = Bridge(T(0.01))
    @test TB.parameters(f) == (; δ)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/566
    outer = Ring(point(6, 4), point(6, 7), point(1, 6), point(1, 1), point(5, 2))
    inner₁ = Ring(point(3, 3), point(3, 4), point(4, 3))
    inner₂ = Ring(point(2, 5), point(2, 6), point(3, 5))
    poly = PolyArea([outer, inner₁, inner₂])
    bpoly = poly |> Bridge(T(0.1))
    @test !hasholes(bpoly)
    @test nvertices(bpoly) == 15

    # unique and bridges
    poly = PolyArea(point.([(0, 0), (1, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1), (0, 1)]))
    cpoly = poly |> Repair{0}() |> Bridge()
    @test cpoly == PolyArea(point.([(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)]))

    # basic ngon tests
    t = Triangle(point(0, 0), point(1, 0), point(0, 1))
    @test (t |> Bridge() |> boundary) == boundary(t)
    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    @test (q |> Bridge() |> boundary) == boundary(q)

    # bridges between holes
    outer = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test vertices(poly) ==
          point.([
      (0, 0),
      (1, 0),
      (1, 1),
      (0, 1),
      (0.2, 0.2),
      (0.2, 0.4),
      (0.4, 0.4),
      (0.4, 0.2),
      (0.6, 0.2),
      (0.6, 0.4),
      (0.8, 0.4),
      (0.8, 0.2)
    ])
    bpoly = poly |> Bridge(T(0.01))
    target =
      point.([
        (-0.0035355339059327372, 0.0035355339059327372),
        (0.19646446609406729, 0.20353553390593274),
        (0.2, 0.4),
        (0.4, 0.405),
        (0.6, 0.405),
        (0.8, 0.4),
        (0.8, 0.2),
        (0.6, 0.2),
        (0.6, 0.395),
        (0.4, 0.395),
        (0.4, 0.2),
        (0.20353553390593274, 0.19646446609406729),
        (0.0035355339059327372, -0.0035355339059327372),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0)
      ])
    @test all(vertices(bpoly) .≈ target)

    poly = Quadrangle(point(0, 1, 0), point(1, 1, 0), point(1, 0, 0), point(0, 0, 0))
    bpoly = poly |> Bridge()
    @test bpoly isa Quadrangle
    @test bpoly == poly
  end

  @testset "Smoothing" begin
    @test !isaffine(LambdaMuSmoothing)
    n, λ, μ = 30, T(0.5), T(0)
    f = LambdaMuSmoothing(n, λ, μ)
    @test TB.parameters(f) == (; n, λ, μ)

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
