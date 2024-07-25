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
    g = cart(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ cart(0, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Rotate(Angle2d(T(π / 2)))
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(0, 0), cart(0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Rotate(Angle2d(T(π / 2)))
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Quadrangle
    @test r ≈ Quadrangle(cart(0, 0), cart(0, 1), cart(-1, 1), cart(-1, 0))
    q = TB.revert(f, r, c)
    @test q isa Quadrangle
    @test q ≈ convert(Quadrangle, g)

    f = Rotate(vector(1, 0, 0), vector(0, 1, 0))
    g = Box(cart(0, 0, 0), cart(1, 1, 1))
    r, c = TB.apply(f, g)
    @test r isa Hexahedron
    @test r ≈ Hexahedron(
      cart(0, 0, 0),
      cart(0, 1, 0),
      cart(-1, 1, 0),
      cart(-1, 0, 0),
      cart(0, 0, 1),
      cart(0, 1, 1),
      cart(-1, 1, 1),
      cart(-1, 0, 1)
    )
    h = TB.revert(f, r, c)
    @test h isa Hexahedron
    @test h ≈ convert(Hexahedron, g)

    # ----------
    # ROPE/RING
    # ----------

    f = Rotate(Angle2d(T(π / 2)))
    g = Rope(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Rope(cart(0, 0), cart(0, 1), cart(-1, 1), cart(-1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Rotate(Angle2d(T(π / 2)))
    g = Ring(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Ring(cart(0, 0), cart(0, 1), cart(-1, 1), cart(-1, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Rotate(AngleAxis(T(π / 2), T(0), T(0), T(1)))
    g = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(cart(0, 0, 0), cart(0, 1, 0), cart(-1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Rotate(vector(0, 0, 1), vector(1, 0, 0))
    g = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(cart(0, 0, 0), cart(0, 0, -1), cart(0, 1, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POLYAREA
    # ---------

    f = Rotate(Angle2d(T(π / 2)))
    p = PolyArea(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    r, c = TB.apply(f, p)
    @test r ≈ PolyArea(cart(0, 0), cart(0, 1), cart(-1, 1), cart(-1, 0))
    @test TB.revert(f, r, c) ≈ p

    # ----------
    # MULTIGEOM
    # ----------

    f = Rotate(Angle2d(T(π / 2)))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Rotate(vector(0, 0, 1), vector(1, 0, 0))
    g = Plane(cart(0, 0, 0), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(cart(0, 0, 0), vector(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Rotate(vector(0, 0, 1), vector(1, 0, 0))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(cart(0, 0, 0), cart(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Rotate(Angle2d(T(π / 2)))
    d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([cart(0, 0), cart(0, 1), cart(-1, 1)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Rotate(Angle2d(T(π / 2)))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
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
    p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
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

    # CRS propagation
    f = Rotate(Angle2d(T(π / 2)))
    p = merc(1, 0)
    @test crs(f(p)) === crs(p)
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
    g = cart(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ cart(2, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Translate(T(1), T(1))
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(1, 1), cart(2, 1))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Translate(T(1), T(1))
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(cart(1, 1), cart(2, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Translate(T(1), T(2), T(3))
    g = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(cart(1, 2, 3), cart(2, 2, 3), cart(1, 3, 4))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Translate(T(1), T(1))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Translate(T(0), T(0), T(1))
    g = Plane(cart(0, 0, 0), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(cart(0, 0, 1), vector(0, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Translate(T(0), T(0), T(1))
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(cart(0, 0, 1), cart(0, 0, 2))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Translate(T(1), T(1))
    d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([cart(1, 1), cart(2, 1), cart(2, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Translate(T(1), T(1))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
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
    @test r ≈ CartesianGrid(cart(1, 1), cart(11, 11), dims=(10, 10))
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
    p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
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
    g = cart(1, 0)
    r, c = TB.apply(f, g)
    @test r ≈ cart(1, 2)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(1, 1), cart(1, 2))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Quadrangle
    @test r ≈ Quadrangle(cart(1, 1), cart(1, 2), cart(0, 2), cart(0, 1))
    q = TB.revert(f, r, c)
    @test q isa Quadrangle
    @test q ≈ convert(Quadrangle, g)

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[1, 2, 3])
    g = Box(cart(0, 0, 0), cart(1, 1, 1))
    r, c = TB.apply(f, g)
    @test r isa Hexahedron
    @test r ≈ Hexahedron(
      cart(1, 2, 3),
      cart(1, 2, 2),
      cart(1, 3, 2),
      cart(1, 3, 3),
      cart(2, 2, 3),
      cart(2, 2, 2),
      cart(2, 3, 2),
      cart(2, 3, 3)
    )
    h = TB.revert(f, r, c)
    @test h isa Hexahedron
    @test h ≈ convert(Hexahedron, g)

    # ---------
    # TRIANGLE
    # ---------

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[1, 2, 3])
    g = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(cart(1, 2, 3), cart(1, 2, 2), cart(2, 3, 3))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[0, 0, 1])
    g = Plane(cart(0, 0, 0), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(cart(0, 0, 1), vector(1, 0, 0))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # CYLINDER
    # ---------

    f = Affine(rotation_between(SVector{3,T}(0, 0, 1), SVector{3,T}(1, 0, 0)), T[0, 0, 1])
    g = Cylinder(T(1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(cart(0, 0, 1), cart(1, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([cart(1, 1), cart(1, 2), cart(0, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
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
    g1 = cart(1, 0)
    g2 = Segment(cart(0, 0), cart(1, 0))
    g3 = Box(cart(0, 0), cart(1, 1))
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

    # CRS propagation
    f = Affine(Angle2d(T(π / 2)), T[1, 1])
    p = merc(1, 0)
    @test crs(f(p)) === crs(p)

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
    g = cart(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ cart(1, 2)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Scale(T(1), T(2))
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(0, 0), cart(1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Scale(T(2), T(1))
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(0, 0), cart(2, 0))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Scale(T(1), T(2))
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(cart(0, 0), cart(1, 2))
    @test TB.revert(f, r, c) ≈ g

    # -----
    # BALL
    # -----

    f = Scale(T(1), T(2))
    g = Ball(cart(1, 2), T(3))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    # -------
    # SPHERE
    # -------

    f = Scale(T(1), T(2))
    g = Sphere(cart(1, 2), T(3))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    f = Scale(T(1), T(2), T(3))
    g = Sphere(cart(1, 2, 3), T(4))
    r, c = TB.apply(f, g)
    @test r isa Ellipsoid
    @test r ≈ Ellipsoid(T.((4, 8, 12)), cart(1, 4, 9))
    m = TB.revert(f, r, c)
    @test m isa SimpleMesh
    @test m ≈ discretize(g)

    f = Scale(T(2))
    g = Sphere(cart(1, 2), T(3))
    r, c = TB.apply(f, g)
    @test r isa Sphere
    @test r ≈ Sphere(cart(2, 4), T(6))
    @test TB.revert(f, r, c) ≈ g

    f = Scale(T(2))
    g = Sphere(cart(1, 2, 3), T(4))
    r, c = TB.apply(f, g)
    @test r isa Sphere
    @test r ≈ Sphere(cart(2, 4, 6), T(8))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # ELLIPSOID
    # ----------

    f = Scale(T(1), T(2), T(3))
    g = Ellipsoid(T.((1, 2, 3)))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    # -----
    # DISK
    # -----

    f = Scale(T(1), T(2), T(3))
    g = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    # -------
    # CIRCLE
    # -------

    f = Scale(T(1), T(2), T(3))
    g = Circle(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    # ----------------
    # CYLINDERSURFACE
    # ----------------

    f = Scale(T(1), T(2), T(3))
    g = CylinderSurface(T(1))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    # ------------------
    # PARABOLOIDSURFACE
    # ------------------

    f = Scale(T(1), T(2), T(3))
    g = ParaboloidSurface(cart(0, 0, 0), T(1), T(2))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test TB.revert(f, r, c) ≈ m

    # ------
    # TORUS
    # ------

    f = Scale(T(1), T(2), T(3))
    g = Torus(cart(1, 1, 1), vector(1, 0, 0), T(2), T(1))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)
    @test centroid(r) ≈ f(centroid(g))
    @test TB.revert(f, r, c) ≈ m

    # ---------
    # TRIANGLE
    # ---------

    f = Scale(T(1), T(2), T(3))
    g = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 2, 3))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Scale(T(1), T(2))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Scale(T(1), T(1), T(2))
    g = Plane(cart(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(cart(1, 1, 2), vector(0, 0, 1))
    @test TB.revert(f, r, c) ≈ g

    f = Scale(T(2), T(1), T(1))
    g = Plane(cart(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Scale(T(1), T(2))
    d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([cart(0, 0), cart(1, 0), cart(1, 2)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Scale(T(1), T(2))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
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
    d = CartesianGrid(cart(1, 1), cart(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(cart(1, 2), cart(11, 22), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Scale(T(1), T(2))
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa RectilinearGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Scale(T(1), T(2))
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa StructuredGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Scale(T(1), T(2))
    p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
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
    g = cart(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ cart(1, 1)
    @test TB.revert(f, r, c) ≈ g

    # --------
    # SEGMENT
    # --------

    f = Stretch(T(1), T(2))
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(0, 0), cart(1, 0))
    @test TB.revert(f, r, c) ≈ g

    f = Stretch(T(2), T(1))
    g = Segment(cart(0, 0), cart(1, 0))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(cart(-0.5, 0), cart(1.5, 0))
    @test TB.revert(f, r, c) ≈ g

    # ----
    # BOX
    # ----

    f = Stretch(T(1), T(2))
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(cart(0, -0.5), cart(1, 1.5))
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # TRIANGLE
    # ---------

    f = Stretch(T(1), T(2), T(2))
    g = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(cart(0, -1 / 3, -1 / 3), cart(1, -1 / 3, -1 / 3), cart(0, 10 / 6, 10 / 6))
    @test TB.revert(f, r, c) ≈ g

    # ----------
    # MULTIGEOM
    # ----------

    f = Stretch(T(1), T(2))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])
    @test TB.revert(f, r, c) ≈ g

    # ------
    # PLANE
    # ------

    f = Stretch(T(1), T(1), T(2))
    g = Plane(cart(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    f = Stretch(T(2), T(1), T(1))
    g = Plane(cart(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ g
    @test TB.revert(f, r, c) ≈ g

    # ---------
    # POINTSET
    # ---------

    f = Stretch(T(1), T(2))
    d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([cart(0, -1 / 3), cart(1, -1 / 3), cart(1, 10 / 6)])
    @test TB.revert(f, r, c) ≈ d

    # ------------
    # GEOMETRYSET
    # ------------

    f = Stretch(T(1), T(2))
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
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
    d = CartesianGrid(cart(1, 1), cart(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(cart(1, -4), cart(11, 16), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Stretch(T(1), T(2))
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa RectilinearGrid
    @test r ≈ SimpleMesh(f(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Stretch(T(1), T(2))
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa StructuredGrid
    @test r ≈ SimpleMesh(f(vertices(d)), topology(d))
    @test TB.revert(f, r, c) ≈ d

    # -----------
    # SIMPLEMESH
    # -----------

    f = Stretch(T(1), T(2))
    p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
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
    d = CartesianGrid(cart(1, 1), cart(11, 11), dims=(10, 10))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid(cart(-0.5, -0.5), cart(0.5, 0.5), dims=(10, 10))
    @test TB.revert(f, r, c) ≈ d

    f = StdCoords()
    d = cartgrid(10, 20)
    r, c = TB.apply(f, d)
    @test r ≈ CartesianGrid(cart(-0.5, -0.5), cart(0.5, 0.5), dims=(10, 20))
    @test TB.revert(f, r, c) ≈ d
    r2 = TB.reapply(f, d, c)
    @test r == r2
  end

  @testset "Proj" begin
    @test !isaffine(Proj(Polar))
    @test !TB.isrevertible(Proj(Polar))
    @test !TB.isinvertible(Proj(Polar))
    @test TB.parameters(Proj(Polar)) == (; CRS=Polar)
    @test TB.parameters(Proj(EPSG{3395})) == (; CRS=Mercator{WGS84Latest})
    @test TB.parameters(Proj(ESRI{54017})) == (; CRS=Behrmann{WGS84Latest})

    # ----
    # VEC
    # ----

    f = Proj(Polar)
    v = vector(1, 0)
    r, c = TB.apply(f, v)
    @test r == v

    # ------
    # POINT
    # ------

    f = Proj(Polar)
    g = cart(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ Point(Polar(T(√2), T(π / 4)))

    # --------
    # SEGMENT
    # --------

    f = Proj(Polar)
    g = Segment(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(Point(Polar(T(0), T(0))), Point(Polar(T(√2), T(π / 4))))

    # ----
    # BOX
    # ----

    f = Proj(Polar)
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(Point(Polar(T(0), T(0))), Point(Polar(T(√2), T(π / 4))))

    # ---------
    # TRIANGLE
    # ---------

    f = Proj(Polar)
    g = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(Point(Polar(T(0), T(0))), Point(Polar(T(1), T(0))), Point(Polar(T(√2), T(π / 4))))

    # ----------
    # MULTIGEOM
    # ----------

    f = Proj(Polar)
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])

    # ------
    # PLANE
    # ------

    f = Proj(Cylindrical)
    g = Plane(cart(1, 1, 1), vector(0, 0, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Plane(Point(Cylindrical(T(√2), T(π / 4), T(1))), vector(0, 0, 1))

    # ---------
    # CYLINDER
    # ---------

    f = Proj(Cylindrical)
    g = Cylinder(cart(0, 0, 0), cart(1, 1, 1))
    r, c = TB.apply(f, g)
    @test r ≈ Cylinder(Point(Cylindrical(T(0), T(0), T(0))), Point(Cylindrical(T(√2), T(π / 4), T(1))))

    # ---------
    # POINTSET
    # ---------

    f = Proj(Polar)
    d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([Point(Polar(T(0), T(0))), Point(Polar(T(1), T(0))), Point(Polar(T(√2), T(π / 4)))])

    # ------------
    # GEOMETRYSET
    # ------------

    f = Proj(Polar)
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])

    # --------------
    # CARTESIANGRID
    # --------------

    f = Proj(Polar)
    d = CartesianGrid((10, 10), cart(1, 1), T.((1, 1)))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid((10, 10), Point(Polar(T(√2), T(π / 4))), T.((1, 1)))

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Proj(Polar)
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa SimpleMesh
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Proj(Polar)
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa SimpleMesh
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

    # -----------
    # SIMPLEMESH
    # -----------

    f = Proj(Polar)
    p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  end

  @testset "LengthUnit" begin
    @test !isaffine(LengthUnit(u"km"))
    @test !TB.isrevertible(LengthUnit(u"cm"))
    @test !TB.isinvertible(LengthUnit(u"km"))
    @test TB.parameters(LengthUnit(u"cm")) == (; unit=u"cm")

    # ----
    # VEC
    # ----

    f = LengthUnit(u"km")
    v = vector(1000, 0)
    r, c = TB.apply(f, v)
    @test r ≈ Vec(T(1) * u"km", T(0) * u"km")

    # ------
    # POINT
    # ------

    f = LengthUnit(u"cm")
    g = cart(1, 1)
    r, c = TB.apply(f, g)
    @test r ≈ Point(T(100) * u"cm", T(100) * u"cm")

    f = LengthUnit(u"km")
    g = Point(Polar(T(1000), T(π / 4)))
    r, c = TB.apply(f, g)
    @test r ≈ Point(Polar(T(1) * u"km", T(π / 4) * u"rad"))

    f = LengthUnit(u"cm")
    g = Point(Cylindrical(T(1), T(π / 4), T(1)))
    r, c = TB.apply(f, g)
    @test r ≈ Point(Cylindrical(T(100) * u"cm", T(π / 4) * u"rad", T(100) * u"cm"))

    f = LengthUnit(u"km")
    g = Point(Spherical(T(1000), T(π / 4), T(π / 4)))
    r, c = TB.apply(f, g)
    @test r ≈ Point(Spherical(T(1) * u"km", T(π / 4) * u"rad", T(π / 4) * u"rad"))

    f = LengthUnit(u"cm")
    g = Point(Mercator(T(1), T(1)))
    @test_throws ArgumentError TB.apply(f, g)

    # --------
    # SEGMENT
    # --------

    f = LengthUnit(u"km")
    g = Segment(cart(0, 0), cart(1000, 1000))
    r, c = TB.apply(f, g)
    @test r ≈ Segment(Point(T(0) * u"km", T(0) * u"km"), Point(T(1) * u"km", T(1) * u"km"))

    # ----
    # BOX
    # ----

    f = LengthUnit(u"cm")
    g = Box(cart(0, 0), cart(1, 1))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r ≈ Box(Point(T(0) * u"cm", T(0) * u"cm"), Point(T(100) * u"cm", T(100) * u"cm"))

    # -------
    # SPHERE
    # -------

    f = LengthUnit(u"km")
    g = Sphere(cart(0, 0), T(1000))
    r, c = TB.apply(f, g)
    @test r isa Sphere
    @test r ≈ Sphere(Point(T(0) * u"km", T(0) * u"km"), T(1) * u"km")

    # ----------
    # ELLIPSOID
    # ----------

    f = LengthUnit(u"cm")
    g = Ellipsoid(T.((1, 1, 1)), cart(0, 0, 0))
    r, c = TB.apply(f, g)
    @test r isa Ellipsoid
    @test r ≈
          Ellipsoid((T(100) * u"cm", T(100) * u"cm", T(100) * u"cm"), Point(T(0) * u"cm", T(0) * u"cm", T(0) * u"cm"))

    # ---------
    # TRIANGLE
    # ---------

    f = LengthUnit(u"km")
    g = Triangle(cart(0, 0), cart(1000, 0), cart(1000, 1000))
    r, c = TB.apply(f, g)
    @test r ≈ Triangle(
      Point(T(0) * u"km", T(0) * u"km"),
      Point(T(1) * u"km", T(0) * u"km"),
      Point(T(1) * u"km", T(1) * u"km")
    )

    # ----------
    # MULTIGEOM
    # ----------

    f = LengthUnit(u"cm")
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r ≈ Multi([f(t), f(t)])

    # ---------
    # POINTSET
    # ---------

    f = LengthUnit(u"km")
    d = PointSet([cart(0, 0), cart(1000, 0), cart(1000, 1000)])
    r, c = TB.apply(f, d)
    @test r ≈ PointSet([
      Point(T(0) * u"km", T(0) * u"km"),
      Point(T(1) * u"km", T(0) * u"km"),
      Point(T(1) * u"km", T(1) * u"km")
    ])

    # ------------
    # GEOMETRYSET
    # ------------

    f = LengthUnit(u"cm")
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r ≈ GeometrySet([f(t), f(t)])
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .≈ [f(t), f(t)])

    # --------------
    # CARTESIANGRID
    # --------------

    f = LengthUnit(u"km")
    d = CartesianGrid((10, 10), cart(1000, 1000), T.((1000, 1000)))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r ≈ CartesianGrid((10, 10), Point(T(1) * u"km", T(1) * u"km"), (T(1) * u"km", T(1) * u"km"))

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = LengthUnit(u"cm")
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa RectilinearGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = LengthUnit(u"km")
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r isa StructuredGrid
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

    # -----------
    # SIMPLEMESH
    # -----------

    f = LengthUnit(u"cm")
    p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  end

  @testset "Shadow" begin
    @test !isaffine(Shadow(:xy))
    @test !TB.isrevertible(Shadow("xy"))
    @test !TB.isinvertible(Shadow(:xy))
    @test TB.parameters(Shadow("xy")) == (; dims=(1, 2))
    @test TB.parameters(Shadow(:yx)) == (; dims=(2, 1))
    @test TB.parameters(Shadow("xz")) == (; dims=(1, 3))
    @test TB.parameters(Shadow(:yz)) == (; dims=(2, 3))
    @test_throws ArgumentError Shadow(:xk)

    # ----
    # VEC
    # ----

    f = Shadow(:xy)
    v = vector(1, 2, 3)
    r, c = TB.apply(f, v)
    @test r == vector(1, 2)

    # ------
    # POINT
    # ------

    f = Shadow(:xz)
    g = cart(1, 2, 3)
    r, c = TB.apply(f, g)
    @test r == cart(1, 3)

    # --------
    # SEGMENT
    # --------

    f = Shadow(:yz)
    g = Segment(cart(1, 2, 3), cart(4, 5, 6))
    r, c = TB.apply(f, g)
    @test r == Segment(cart(2, 3), cart(5, 6))

    # ----
    # BOX
    # ----

    f = Shadow(:xy)
    g = Box(cart(1, 2, 3), cart(4, 5, 6))
    r, c = TB.apply(f, g)
    @test r isa Box
    @test r == Box(cart(1, 2), cart(4, 5))

    # ------
    # PLANE
    # ------

    f = Shadow(:xz)
    g = Plane(cart(0, 0, 0), vector(0, 0, 1))
    @test_throws ArgumentError TB.apply(f, g)

    # ----------
    # ELLIPSOID
    # ----------

    f = Shadow(:yz)
    g = Ellipsoid(T.((1, 2, 3)))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # -----
    # DISK
    # -----

    f = Shadow(:xy)
    g = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # -------
    # CIRCLE
    # -------

    f = Shadow(:xz)
    g = Circle(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # ----------------
    # CYLINDERSURFACE
    # ----------------

    f = Shadow(:yz)
    g = CylinderSurface(T(1))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # ------------
    # CONESURFACE
    # ------------

    f = Shadow(:xy)
    p = Plane(cart(0, 0, 0), vector(0, 0, 1))
    g = ConeSurface(Disk(p, T(2)), cart(0, 0, 1))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # ---------------
    # FRUSTUMSURFACE
    # ---------------

    f = Shadow(:xz)
    pb = Plane(cart(0, 0, 0), vector(0, 0, 1))
    pt = Plane(cart(0, 0, 10), vector(0, 0, 1))
    g = FrustumSurface(Disk(pb, T(1)), Disk(pt, T(2)))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # ------------------
    # PARABOLOIDSURFACE
    # ------------------

    f = Shadow(:yz)
    g = ParaboloidSurface(cart(0, 0, 0), T(1), T(2))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # ------
    # TORUS
    # ------

    f = Shadow(:xy)
    g = Torus(cart(1, 1, 1), vector(1, 0, 0), T(2), T(1))
    m = discretize(g)
    r, c = TB.apply(f, g)
    @test r isa SimpleMesh
    @test r ≈ f(m)

    # ---------
    # TRIANGLE
    # ---------

    f = Shadow(:xz)
    g = Triangle(cart(1, 2, 3), cart(4, 5, 6), cart(7, 8, 9))
    r, c = TB.apply(f, g)
    @test r == Triangle(cart(1, 3), cart(4, 6), cart(7, 9))

    # ----------
    # MULTIGEOM
    # ----------

    f = Shadow(:yz)
    t = Triangle(cart(1, 2, 3), cart(4, 5, 6), cart(7, 8, 9))
    g = Multi([t, t])
    r, c = TB.apply(f, g)
    @test r == Multi([f(t), f(t)])

    # ---------
    # POINTSET
    # ---------

    f = Shadow(:xy)
    d = PointSet([cart(1, 2, 3), cart(4, 5, 6), cart(7, 8, 9)])
    r, c = TB.apply(f, d)
    @test r == PointSet([cart(1, 2), cart(4, 5), cart(7, 8)])

    # ------------
    # GEOMETRYSET
    # ------------

    f = Shadow(:xz)
    t = Triangle(cart(1, 2, 3), cart(4, 5, 6), cart(7, 8, 9))
    d = GeometrySet([t, t])
    r, c = TB.apply(f, d)
    @test r == GeometrySet([f(t), f(t)])
    d = [t, t]
    r, c = TB.apply(f, d)
    @test all(r .== [f(t), f(t)])

    # --------------
    # CARTESIANGRID
    # --------------

    f = Shadow(:yz)
    d = CartesianGrid((10, 11, 12), cart(1, 2, 3), T.((1.0, 1.1, 1.2)))
    r, c = TB.apply(f, d)
    @test r isa CartesianGrid
    @test r == CartesianGrid((11, 12), cart(2, 3), T.((1.1, 1.2)))

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Shadow(:xy)
    d = convert(RectilinearGrid, cartgrid(10, 11, 12))
    r, c = TB.apply(f, d)
    @test r isa RectilinearGrid
    @test r == RectilinearGrid(Meshes.xyz(d)[1], Meshes.xyz(d)[2])

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Shadow(:xz)
    d = convert(StructuredGrid, cartgrid(10, 11, 12))
    r, c = TB.apply(f, d)
    @test r isa StructuredGrid
    @test r == StructuredGrid(Meshes.XYZ(d)[1][:, 1, :], Meshes.XYZ(d)[3][:, 1, :])

    # -----------
    # SIMPLEMESH
    # -----------

    f = Shadow(:yz)
    p = cart.([(0, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (0, 0.5, 0.5)])
    c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    d = SimpleMesh(p, c)
    r, c = TB.apply(f, d)
    @test r == SimpleMesh(f.(vertices(d)), topology(d))
  end

  @testset "Within" begin
    @test !isaffine(Within(x=(T(2), T(4))))
    @test !TB.isrevertible(Within(x=(T(2), T(4))))
    @test !TB.isinvertible(Within(x=(T(2), T(4))))
    @test TB.parameters(Within(x=(T(2), T(4)))) == (x=(T(2) * u"m", T(4) * u"m"), y=nothing, z=nothing)
    @test TB.parameters(Within(y=(T(2) * u"km", T(4) * u"km"))) ==
          (x=nothing, y=(T(2) * u"km", T(4) * u"km"), z=nothing)
    @test TB.parameters(Within(z=(2, 4))) == (x=nothing, y=nothing, z=(2.0u"m", 4.0u"m"))
    @test_throws ArgumentError Within(x=(T(2) * u"°", T(4) * u"°"))

    # ---------
    # POINTSET
    # ---------

    f = Within(x=(T(1.5), T(3.5)))
    d = PointSet([cart(1, 0), cart(2, 1), cart(3, 1), cart(4, 0)])
    r, c = TB.apply(f, d)
    @test r == PointSet([cart(2, 1), cart(3, 1)])

    # ------------
    # GEOMETRYSET
    # ------------

    f = Within(x=(T(1.5), T(3.5)))
    t1 = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    t2 = t1 |> Translate(T(2), T(2))
    t3 = t2 |> Translate(T(2), T(2))
    d = GeometrySet([t1, t2, t3])
    r, c = TB.apply(f, d)
    @test r == GeometrySet([t2])

    # --------------
    # CARTESIANGRID
    # --------------

    f = Within(z=(T(1.5), T(4.5)))
    d = cartgrid(10, 10, 10)
    r, c = TB.apply(f, d)
    @test r == CartesianGrid((10, 10, 4), cart(0, 0, 1), T.((1, 1, 1)))

    # ----------------
    # RECTILINEARGRID
    # ----------------

    f = Within(y=(T(3.5), T(6.5)))
    d = convert(RectilinearGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r == convert(RectilinearGrid, CartesianGrid((10, 4), cart(0, 3), T.((1, 1))))

    # ---------------
    # STRUCTUREDGRID
    # ---------------

    f = Within(x=(T(5.5), T(8.5)))
    d = convert(StructuredGrid, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r == convert(StructuredGrid, CartesianGrid((4, 10), cart(5, 0), T.((1, 1))))

    # -----------
    # SIMPLEMESH
    # -----------

    f = Within(x=(T(1.5), T(4.5)), y=(T(3.5), T(6.5)))
    d = convert(SimpleMesh, cartgrid(10, 10))
    r, c = TB.apply(f, d)
    @test r == convert(SimpleMesh, CartesianGrid((4, 4), cart(1, 3), T.((1, 1))))
  end

  @testset "Repair{0}" begin
    @test !isaffine(Repair)
    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 0), (1, 1), (0, 1), (0, 1)]))
    rpoly = poly |> Repair{0}()
    @test nvertices(rpoly) == 4
    @test vertices(rpoly) == cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
  end

  @testset "Repair{1}" begin
    # a tetrahedron with an unused vertex
    points = cart.([(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)])
    connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)])
    mesh = SimpleMesh(points, connec)
    rmesh = mesh |> Repair{1}()
    @test nvertices(rmesh) == nvertices(mesh) - 1
    @test nelements(rmesh) == nelements(mesh)
    @test cart(5, 5, 5) ∉ vertices(rmesh)
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
      cart.([(0.0, 0.0), (0.5, -0.5), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (-0.5, 0.5)])
    )
    rpoly = poly |> Repair{8}()
    @test nvertices(rpoly) == 4
    @test vertices(rpoly) == cart.([(0.5, -0.5), (1.5, 0.5), (0.5, 1.5), (-0.5, 0.5)])

    # degenerate triangle with repeated vertices
    poly = PolyArea(cart.([(0, 0), (1, 1), (1, 1)]))
    rpoly = poly |> Repair{8}()
    @test !hasholes(rpoly)
    @test rings(rpoly) == [Ring(cart(0, 0))]
    @test vertices(rpoly) == [cart(0, 0)]
  end

  @testset "Repair{9}" begin
    poly = Quadrangle(cart(0, 1, 0), cart(1, 1, 0), cart(1, 0, 0), cart(0, 0, 0))
    bpoly = poly |> Repair{9}()
    @test bpoly isa Quadrangle
    @test bpoly == poly
  end

  @testset "Repair{10}" begin
    outer = Ring(cart.([(0, 0), (0, 3), (2, 3), (2, 2), (3, 2), (3, 0)]))
    inner = Ring(cart.([(1, 1), (1, 2), (2, 2), (2, 1)]))
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
    outer = Ring(cart(6, 4), cart(6, 7), cart(1, 6), cart(1, 1), cart(5, 2))
    inner₁ = Ring(cart(3, 3), cart(3, 4), cart(4, 3))
    inner₂ = Ring(cart(2, 5), cart(2, 6), cart(3, 5))
    poly = PolyArea([outer, inner₁, inner₂])
    bpoly = poly |> Bridge(T(0.1))
    @test !hasholes(bpoly)
    @test nvertices(bpoly) == 15

    # unique and bridges
    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1), (0, 1)]))
    cpoly = poly |> Repair{0}() |> Bridge()
    @test cpoly == PolyArea(cart.([(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)]))

    # basic ngon tests
    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test (t |> Bridge() |> boundary) == boundary(t)
    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test (q |> Bridge() |> boundary) == boundary(q)

    # bridges between holes
    outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test vertices(poly) ==
          cart.([
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
      cart.([
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

    poly = Quadrangle(cart(0, 1, 0), cart(1, 1, 0), cart(1, 0, 0), cart(0, 0, 0))
    bpoly = poly |> Bridge()
    @test bpoly isa Quadrangle
    @test bpoly == poly

    # bridge with latlon coords
    outer = latlon.([(0, 0), (0, 90), (90, 90), (90, 0)])
    hole1 = latlon.([(10, 10), (10, 20), (20, 20), (20, 10)])
    hole2 = latlon.([(10, 80), (10, 90), (20, 90), (20, 80)])
    poly = PolyArea([outer, hole1, hole2])
    bpoly = poly |> Bridge()
    @test nvertices(bpoly) == 16
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
