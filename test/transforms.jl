@testitem "Rotate" setup = [Setup] begin
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
  @test r ≈ Quadrangle(cart(0, 0), cart(0, 1), cart(-1, 1), cart(-1, 0))
  @test TB.revert(f, r, c) ≈ g

  f = Rotate(vector(1, 0, 0), vector(0, 1, 0))
  g = Box(cart(0, 0, 0), cart(1, 1, 1))
  r, c = TB.apply(f, g)
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

  # ----------
  # ELLIPSOID
  # ----------

  R = RotXYZ(T(π / 4), T(π / 5), T(π / 3))
  f = Rotate(R)
  g = Ellipsoid((T(3), T(2), T(1)), (T(1), T(1), T(1)))
  r, c = TB.apply(f, g)
  @test center(r) == center(g) |> Rotate(R)
  @test rotation(r) == R

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
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  @test TB.revert(f, r, c) ≈ d

  # ------------
  # REGULARGRID
  # ------------

  f = Rotate(Angle2d(T(π / 2)))
  d = RegularGrid(merc(0, 0), merc(1, 1), dims=(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  @test TB.revert(f, r, c) ≈ d

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = Rotate(Angle2d(T(π / 2)))
  d = convert(RectilinearGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  @test TB.revert(f, r, c) ≈ d

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = Rotate(Angle2d(T(π / 2)))
  d = convert(StructuredGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
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

@testitem "Translate" setup = [Setup] begin
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

  # ------------
  # REGULARGRID
  # ------------

  f = Translate(T(1), T(1))
  d = RegularGrid((8, 8), Point(Polar(T(0), T(0))), (T(1), T(π / 4)))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  @test TB.revert(f, r, c) ≈ d

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

  f = Translate(T(1), T(1))
  g = RegularGrid((8, 8), Point(Polar(T(0), T(0))), (T(1), T(π / 4)))
  d = convert(RectilinearGrid, g)
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
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

  f = Translate(T(1), T(1))
  g = RegularGrid((8, 8), Point(Polar(T(0), T(0))), (T(1), T(π / 4)))
  d = convert(StructuredGrid, g)
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
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

@testitem "Affine" setup = [Setup] begin
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
  @test r ≈ Quadrangle(cart(1, 1), cart(1, 2), cart(0, 2), cart(0, 1))
  @test TB.revert(f, r, c) ≈ g

  f = Affine(Meshes.urotbetween(vector(0, 0, 1), vector(1, 0, 0)), T[1, 2, 3])
  g = Box(cart(0, 0, 0), cart(1, 1, 1))
  r, c = TB.apply(f, g)
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
  @test h ≈ convert(Hexahedron, g)

  # -----
  # BALL
  # -----

  f = Affine(Diagonal(T[1, 2]), T[1, 1])
  g = Ball(cart(1, 2), T(3))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # -------
  # SPHERE
  # -------

  f = Affine(Diagonal(T[1, 2]), T[1, 1])
  g = Sphere(cart(1, 2), T(3))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ----------
  # ELLIPSOID
  # ----------

  f = Affine(Diagonal(T[1, 2, 3]), T[1, 1, 1])
  g = Ellipsoid(T.((4, 5, 6)), cart(1, 2, 3))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # -----
  # DISK
  # -----

  f = Affine(Diagonal(T[1, 2, 3]), T[1, 1, 1])
  g = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # -------
  # CIRCLE
  # -------

  f = Affine(Diagonal(T[1, 2, 3]), T[1, 1, 1])
  g = Circle(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ----------------
  # CYLINDERSURFACE
  # ----------------

  f = Affine(Diagonal(T[1, 2, 3]), T[1, 1, 1])
  g = CylinderSurface(T(1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ------------------
  # PARABOLOIDSURFACE
  # ------------------

  f = Affine(Diagonal(T[1, 2, 3]), T[1, 1, 1])
  g = ParaboloidSurface(cart(0, 0, 0), T(1), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ------
  # TORUS
  # ------

  f = Affine(Diagonal(T[1, 2, 3]), T[1, 1, 1])
  g = Torus(cart(1, 1, 1), vector(1, 0, 0), T(2), T(1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ---------
  # TRIANGLE
  # ---------

  f = Affine(Meshes.urotbetween(vector(0, 0, 1), vector(1, 0, 0)), T[1, 2, 3])
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

  f = Affine(Meshes.urotbetween(vector(0, 0, 1), vector(1, 0, 0)), T[0, 0, 1])
  g = Plane(cart(0, 0, 0), vector(0, 0, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Plane(cart(0, 0, 1), vector(1, 0, 0))
  @test TB.revert(f, r, c) ≈ g

  # ---------
  # CYLINDER
  # ---------

  f = Affine(Meshes.urotbetween(vector(0, 0, 1), vector(1, 0, 0)), T[0, 0, 1])
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
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  @test TB.revert(f, r, c) ≈ d

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = Affine(Angle2d(T(π / 2)), T[1, 1])
  d = convert(RectilinearGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
  @test TB.revert(f, r, c) ≈ d

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = Affine(Angle2d(T(π / 2)), T[1, 1])
  d = convert(StructuredGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
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

@testitem "Scale" setup = [Setup] begin
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
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  f = Scale(T(2))
  g = Ball(cart(1, 2), T(3))
  r, c = TB.apply(f, g)
  @test r isa Ball
  @test r ≈ Ball(cart(2, 4), T(6))
  @test TB.revert(f, r, c) ≈ g

  f = Scale(T(2))
  g = Ball(cart(1, 2, 3), T(4))
  r, c = TB.apply(f, g)
  @test r isa Ball
  @test r ≈ Ball(cart(2, 4, 6), T(8))
  @test TB.revert(f, r, c) ≈ g

  # -------
  # SPHERE
  # -------

  f = Scale(T(1), T(2))
  g = Sphere(cart(1, 2), T(3))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  f = Scale(T(1), T(2), T(3))
  g = Sphere(cart(1, 2, 3), T(4))
  r, c = TB.apply(f, g)
  @test r isa Ellipsoid
  @test r ≈ Ellipsoid(T.((4, 8, 12)), cart(1, 4, 9))
  @test discretize(TB.revert(f, r, c)) ≈ discretize(g)

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
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  f = Scale(T(2))
  g = Ellipsoid(T.((4, 5, 6)), cart(1, 2, 3))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test r isa Ellipsoid
  @test r ≈ Ellipsoid(T.((8, 10, 12)), cart(2, 4, 6))
  @test TB.revert(f, r, c) ≈ g

  # -----
  # DISK
  # -----

  f = Scale(T(1), T(2), T(3))
  g = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # -------
  # CIRCLE
  # -------

  f = Scale(T(1), T(2), T(3))
  g = Circle(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ----------------
  # CYLINDERSURFACE
  # ----------------

  f = Scale(T(1), T(2), T(3))
  g = CylinderSurface(T(1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ------------------
  # PARABOLOIDSURFACE
  # ------------------

  f = Scale(T(1), T(2), T(3))
  g = ParaboloidSurface(cart(0, 0, 0), T(1), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test discretize(TB.revert(f, r, c)) ≈ m

  # ------
  # TORUS
  # ------

  f = Scale(T(1), T(2), T(3))
  g = Torus(cart(1, 1, 1), vector(1, 0, 0), T(2), T(1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)
  @test centroid(r) ≈ f(centroid(g))
  @test discretize(TB.revert(f, r, c)) ≈ m

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

  # ------------
  # REGULARGRID
  # ------------

  f = Scale(T(1), T(2))
  d = RegularGrid(merc(1, 1), merc(11, 11), dims=(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ RegularGrid(merc(1, 2), merc(11, 22), dims=(10, 10))
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

@testitem "Stretch" setup = [Setup] begin
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

@testitem "StdCoords" setup = [Setup] begin
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

@testitem "Proj" setup = [Setup] begin
  @test !isaffine(Proj(Polar))
  @test !TB.isrevertible(Proj(Polar))
  @test !TB.isinvertible(Proj(Polar))
  @test TB.parameters(Proj(Polar)) == (; CRS=Polar)
  @test TB.parameters(Proj(EPSG{3395})) == (; CRS=Mercator{WGS84Latest})
  @test TB.parameters(Proj(ESRI{54017})) == (; CRS=Behrmann{WGS84Latest})
  f = Proj(Mercator)
  @test sprint(show, f) == "Proj(CRS: CoordRefSystems.Mercator)"
  @test sprint(show, MIME"text/plain"(), f) == """
  Proj transform
  └─ CRS: CoordRefSystems.Mercator"""

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

  # change the manifold
  f = Proj(Mercator)
  g = latlon(0, 0)
  r, c = TB.apply(f, g)
  @test manifold(r) === 𝔼{2}
  @test r ≈ merc(0, 0)

  # preserve the manifold
  f = Proj(Cartesian)
  g = latlon(0, 0)
  r, c = TB.apply(f, g)
  @test manifold(r) === 🌐
  @test r ≈ Point(Cartesian{WGS84Latest}(T(6378137.0), T(0), T(0)))

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

  # ---------
  # CYLINDER
  # ---------

  f = Proj(Cylindrical)
  g = Cylinder(cart(0, 0, 0), cart(1, 1, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Cylinder(Point(Cylindrical(T(0), T(0), T(0))), Point(Cylindrical(T(√2), T(π / 4), T(1))))

  # --------------------
  # TRANSFORMEDGEOMETRY
  # --------------------

  f = Proj(Mercator)
  b = Box(latlon(0, 0), latlon(45, 45))
  g = TransformedGeometry(b, Identity())
  r, c = TB.apply(f, g)
  @test r ≈ TransformedGeometry(b, f)

  f = Proj(LatLon)
  b = Box(merc(0, 0), merc(1, 1))
  g = TransformedGeometry(b, Identity())
  r, c = TB.apply(f, g)
  @test r ≈ TransformedGeometry(b, f)

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
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = Proj(Polar)
  d = convert(RectilinearGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = Proj(Polar)
  d = convert(StructuredGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
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

  # ----------
  # SUBDOMAIN
  # ----------

  f = Proj(Polar)
  g = CartesianGrid((10, 10), cart(1, 1), T.((1, 1)))
  d = view(g, 1:10)
  r, c = TB.apply(f, d)
  @test r ≈ view(SimpleMesh(f.(vertices(g)), topology(g)), 1:10)

  # --------------
  # SPECIAL CASES
  # --------------

  f = Proj(Mercator)
  g = Box(latlon(0, 180), latlon(45, 90))
  r, c = TB.apply(f, g)
  @test manifold(r) === 𝔼{2}

  f = Proj(LatLon)
  g = Box(merc(0, 0), merc(1, 1))
  r, c = TB.apply(f, g)
  @test manifold(r) === 🌐

  # --------------
  # NO CONVERSION
  # --------------

  f = Proj(Cartesian)
  g = cart(1, 1)
  r, c = TB.apply(f, g)
  @test r === g
  f = Proj(crs(cart(0, 0)))
  r, c = TB.apply(f, g)
  @test r === g

  f = Proj(LatLon)
  g = Ring(latlon(0, 0), latlon(0, 1), latlon(1, 0))
  r, c = TB.apply(f, g)
  @test r === g
  f = Proj(crs(latlon(0, 0)))
  r, c = TB.apply(f, g)
  @test r === g

  f = Proj(Mercator)
  p = merc.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  d = SimpleMesh(p, c)
  r, c = TB.apply(f, d)
  @test r === d
  f = Proj(crs(merc(0, 0)))
  r, c = TB.apply(f, d)
  @test r === d
end

@testitem "Morphological" setup = [Setup] begin
  f = Morphological(c -> c)
  @test !isaffine(f)
  @test !TB.isrevertible(f)
  @test !TB.isinvertible(f)
  @test TB.parameters(f) == (; fun=f.fun)

  # ----
  # VEC
  # ----

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  v = vector(1, 0)
  r, c = TB.apply(f, v)
  @test r == v

  # ------
  # POINT
  # ------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  g = cart(1, 1)
  r, c = TB.apply(f, g)
  @test r == cart(1, 1, 0)

  # --------
  # SEGMENT
  # --------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  g = Segment(cart(0, 0), cart(1, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Segment(cart(0, 0, 0), cart(1, 1, 0))

  # ----
  # BOX
  # ----

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  g = Box(cart(0, 0), cart(1, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))

  # ---------
  # TRIANGLE
  # ---------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  g = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0))

  # ----------
  # MULTIGEOM
  # ----------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  g = Multi([t, t])
  r, c = TB.apply(f, g)
  @test r ≈ Multi([f(t), f(t)])

  # --------------------
  # TRANSFORMEDGEOMETRY
  # --------------------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  b = Box(cart(0, 0), cart(1, 1))
  g = TransformedGeometry(b, Identity())
  r, c = TB.apply(f, g)
  @test r ≈ Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))

  # ---------
  # POINTSET
  # ---------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
  r, c = TB.apply(f, d)
  @test r ≈ PointSet([cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0)])

  # ------------
  # GEOMETRYSET
  # ------------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
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

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  d = CartesianGrid((10, 10), cart(1, 1), T.((1, 1)))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  d = convert(RectilinearGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  d = convert(StructuredGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # -----------
  # SIMPLEMESH
  # -----------

  f = Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
  p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  d = SimpleMesh(p, c)
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))
end

@testitem "ReinterpretCoords" setup = [Setup] begin
  f = ReinterpretCoords(Cartesian, LatLon)
  @test TB.isrevertible(f)
  @test TB.isinvertible(f)
  @test TB.parameters(f) == (; CRS₁=Cartesian, CRS₂=LatLon)
  @test inverse(f) == ReinterpretCoords(LatLon, Cartesian)

  # ----
  # VEC
  # ----

  f = ReinterpretCoords(Cartesian, LatLon)
  v = vector(1, 0)
  r, c = TB.apply(f, v)
  @test r == v

  # ------
  # POINT
  # ------

  f = ReinterpretCoords(Cartesian, LatLon)
  g = cart(1, 2)
  r, c = TB.apply(f, g)
  @test r == latlon(2, 1)

  # --------
  # SEGMENT
  # --------

  f = ReinterpretCoords(Cartesian, LatLon)
  g = Segment(cart(1, 2), cart(3, 4))
  r, c = TB.apply(f, g)
  @test r ≈ Segment(latlon(2, 1), latlon(4, 3))

  # ----
  # BOX
  # ----

  f = ReinterpretCoords(Cartesian, LatLon)
  g = Box(cart(0, 0), cart(1, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Quadrangle(latlon(0, 0), latlon(0, 1), latlon(1, 1), latlon(1, 0))

  # ---------
  # TRIANGLE
  # ---------

  f = ReinterpretCoords(Cartesian, LatLon)
  g = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  r, c = TB.apply(f, g)
  @test r ≈ Triangle(latlon(0, 0), latlon(0, 1), latlon(1, 1))

  # ----------
  # MULTIGEOM
  # ----------

  f = ReinterpretCoords(Cartesian, LatLon)
  t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  g = Multi([t, t])
  r, c = TB.apply(f, g)
  @test r ≈ Multi([f(t), f(t)])

  # --------------------
  # TRANSFORMEDGEOMETRY
  # --------------------

  f = ReinterpretCoords(Cartesian, LatLon)
  b = Box(cart(0, 0), cart(1, 1))
  g = TransformedGeometry(b, Identity())
  r, c = TB.apply(f, g)
  @test r ≈ Quadrangle(latlon(0, 0), latlon(0, 1), latlon(1, 1), latlon(1, 0))

  # ---------
  # POINTSET
  # ---------

  f = ReinterpretCoords(Cartesian, LatLon)
  d = PointSet([cart(0, 0), cart(1, 0), cart(1, 1)])
  r, c = TB.apply(f, d)
  @test r ≈ PointSet([latlon(0, 0), latlon(0, 1), latlon(1, 1)])

  # ------------
  # GEOMETRYSET
  # ------------

  f = ReinterpretCoords(Cartesian, LatLon)
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

  f = ReinterpretCoords(Cartesian, LatLon)
  d = cartgrid(10, 10)
  r, c = TB.apply(f, d)
  @test r isa Grid
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = ReinterpretCoords(Cartesian, LatLon)
  d = convert(RectilinearGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r isa Grid
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = ReinterpretCoords(Cartesian, LatLon)
  d = convert(StructuredGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r isa Grid
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # -----------
  # SIMPLEMESH
  # -----------

  f = ReinterpretCoords(Cartesian, LatLon)
  p = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
  c = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  d = SimpleMesh(p, c)
  r, c = TB.apply(f, d)
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # throws error with invalid arguments
  f = ReinterpretCoords(Cartesian, LatLon)
  @test_throws ArgumentError f(latlon(0, 0))
  f = ReinterpretCoords(LatLon, Cartesian)
  @test_throws ArgumentError f(cart(0, 0))
end

@testitem "LengthUnit" setup = [Setup] begin
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
  @test r ≈ Ellipsoid((T(100) * u"cm", T(100) * u"cm", T(100) * u"cm"), Point(T(0) * u"cm", T(0) * u"cm", T(0) * u"cm"))

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

  # ------------
  # REGULARGRID
  # ------------

  f = LengthUnit(u"cm")
  d = RegularGrid((8, 8), Point(Polar(T(1), T(0))), (T(1), T(π / 4)))
  r, c = TB.apply(f, d)
  @test r ≈ RegularGrid((8, 8), Point(Polar(T(100) * u"cm", T(0) * u"rad")), (T(100) * u"cm", T(π / 4) * u"rad"))

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
  @test r ≈ SimpleMesh(f.(vertices(d)), topology(d))

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = LengthUnit(u"km")
  d = convert(StructuredGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
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

@testitem "ValidCoords" setup = [Setup] begin
  @test !isaffine(ValidCoords(Mercator))
  @test !TB.isrevertible(ValidCoords(Mercator))
  @test !TB.isinvertible(ValidCoords(Mercator))
  @test TB.parameters(ValidCoords(Mercator)) == (; CRS=Mercator)
  @test TB.parameters(ValidCoords(EPSG{3395})) == (; CRS=Mercator{WGS84Latest})
  @test TB.parameters(ValidCoords(ESRI{54017})) == (; CRS=Behrmann{WGS84Latest})
  f = ValidCoords(Mercator)
  @test sprint(show, f) == "ValidCoords(CRS: CoordRefSystems.Mercator)"
  @test sprint(show, MIME"text/plain"(), f) == """
  ValidCoords transform
  └─ CRS: CoordRefSystems.Mercator"""

  # ---------
  # POINTSET
  # ---------

  f = ValidCoords(Mercator)
  d = PointSet([latlon(-90, 0), latlon(-45, 0), latlon(0, 0), latlon(45, 0), latlon(90, 0)])
  r, c = TB.apply(f, d)
  @test r == PointSet([latlon(-45, 0), latlon(0, 0), latlon(45, 0)])

  # ------------
  # GEOMETRYSET
  # ------------

  f = ValidCoords(Mercator)
  t1 = Triangle(latlon(-90, 0), latlon(-45, 30), latlon(-45, -30))
  t2 = Triangle(latlon(-45, 0), latlon(0, 30), latlon(0, -30))
  t3 = Triangle(latlon(0, -30), latlon(0, 30), latlon(45, 0))
  t4 = Triangle(latlon(45, -30), latlon(45, 30), latlon(90, 0))
  d = GeometrySet([t1, t2, t3, t4])
  r, c = TB.apply(f, d)
  @test r == GeometrySet([t2, t3])

  f = ValidCoords(Mercator)
  t1 = Triangle(latlon(-90, 0), latlon(-45, 30), latlon(-45, -30))
  t2 = Triangle(latlon(-45, 0), latlon(0, 30), latlon(0, -30))
  t3 = Triangle(latlon(0, -30), latlon(0, 30), latlon(45, 0))
  t4 = Triangle(latlon(45, -30), latlon(45, 30), latlon(90, 0))
  m1 = Multi([t1, t4])
  m2 = Multi([t2, t3])
  d = GeometrySet([m1, m2])
  r, c = TB.apply(f, d)
  @test r == GeometrySet([m2])

  f = ValidCoords(Mercator)
  b1 = Box(latlon(-90, -30), latlon(-30, 30))
  b2 = Box(latlon(-30, -30), latlon(30, 30))
  b3 = Box(latlon(30, -30), latlon(90, 30))
  d = GeometrySet([b1, b2, b3])
  r, c = TB.apply(f, d)
  @test r == GeometrySet([b2])

  # --------------
  # CARTESIANGRID
  # --------------

  f = ValidCoords(Mercator)
  d = RegularGrid(latlon(-90, -180), latlon(90, 180), dims=(3, 3))
  r, c = TB.apply(f, d)
  linds = LinearIndices(size(d))
  @test r == view(d, [linds[2, 1], linds[2, 2], linds[2, 3]])

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = ValidCoords(Mercator)
  g = RegularGrid(latlon(-90, -180), latlon(90, 180), dims=(3, 3))
  d = convert(RectilinearGrid, g)
  r, c = TB.apply(f, d)
  linds = LinearIndices(size(d))
  @test r == view(d, [linds[2, 1], linds[2, 2], linds[2, 3]])

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = ValidCoords(Mercator)
  g = RegularGrid(latlon(-90, -180), latlon(90, 180), dims=(3, 3))
  d = convert(StructuredGrid, g)
  r, c = TB.apply(f, d)
  linds = LinearIndices(size(d))
  @test r == view(d, [linds[2, 1], linds[2, 2], linds[2, 3]])

  # -----------
  # SIMPLEMESH
  # -----------

  f = ValidCoords(Mercator)
  g = RegularGrid(latlon(-90, -180), latlon(90, 180), dims=(3, 3))
  d = convert(SimpleMesh, g)
  r, c = TB.apply(f, d)
  linds = LinearIndices(size(d))
  @test r == view(d, [linds[2, 1], linds[2, 2], linds[2, 3]])
end

@testitem "Shadow" setup = [Setup] begin
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
  @test discretize(r) ≈ f(m)

  # -----
  # DISK
  # -----

  f = Shadow(:xy)
  g = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

  # -------
  # CIRCLE
  # -------

  f = Shadow(:xz)
  g = Circle(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

  # ----------------
  # CYLINDERSURFACE
  # ----------------

  f = Shadow(:yz)
  g = CylinderSurface(T(1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

  # ------------
  # CONESURFACE
  # ------------

  f = Shadow(:xy)
  p = Plane(cart(0, 0, 0), vector(0, 0, 1))
  g = ConeSurface(Disk(p, T(2)), cart(0, 0, 1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

  # ---------------
  # FRUSTUMSURFACE
  # ---------------

  f = Shadow(:xz)
  pb = Plane(cart(0, 0, 0), vector(0, 0, 1))
  pt = Plane(cart(0, 0, 10), vector(0, 0, 1))
  g = FrustumSurface(Disk(pb, T(1)), Disk(pt, T(2)))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

  # ------------------
  # PARABOLOIDSURFACE
  # ------------------

  f = Shadow(:yz)
  g = ParaboloidSurface(cart(0, 0, 0), T(1), T(2))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

  # ------
  # TORUS
  # ------

  f = Shadow(:xy)
  g = Torus(cart(1, 1, 1), vector(1, 0, 0), T(2), T(1))
  m = discretize(g)
  r, c = TB.apply(f, g)
  @test discretize(r) ≈ f(m)

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

  # ------------
  # REGULARGRID
  # ------------

  f = Shadow(:yz)
  d = RegularGrid((8, 8, 8), Point(Cylindrical(T(0), T(0), T(0))), (T(1), T(π / 4), T(1)))
  r, c = TB.apply(f, d)
  @test r == SimpleMesh(f.(vertices(d)), topology(d))

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
  @test r == SimpleMesh(f.(vertices(d)), topology(d))

  # ---------------
  # STRUCTUREDGRID
  # ---------------

  f = Shadow(:xz)
  d = convert(StructuredGrid, cartgrid(10, 11, 12))
  r, c = TB.apply(f, d)
  @test r == SimpleMesh(f.(vertices(d)), topology(d))

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

@testitem "Slice" setup = [Setup] begin
  @test !isaffine(Slice(x=(T(2), T(4))))
  @test !TB.isrevertible(Slice(x=(T(2), T(4))))
  @test !TB.isinvertible(Slice(x=(T(2), T(4))))
  @test TB.parameters(Slice(x=(T(2), T(4)))) == (; limits=(; x=(T(2), T(4))))
  @test TB.parameters(Slice(y=(T(2) * u"km", T(4) * u"km"))) == (; limits=(; y=(T(2) * u"km", T(4) * u"km")))
  @test TB.parameters(Slice(z=(2, 4))) == (; limits=(; z=(2, 4)))
  @test TB.parameters(Slice(lat=(30, 60))) == (; limits=(; lat=(30, 60)))
  @test TB.parameters(Slice(lon=(45u"°", 90u"°"))) == (; limits=(; lon=(45u"°", 90u"°")))

  # --------------
  # CARTESIANGRID
  # --------------

  f = Slice(z=(T(1.5), T(4.5)))
  d = cartgrid(10, 10, 10)
  r, c = TB.apply(f, d)
  @test r isa CartesianGrid
  @test r == CartesianGrid((10, 10, 4), cart(0, 0, 1), T.((1, 1, 1)))

  # ----------------
  # RECTILINEARGRID
  # ----------------

  f = Slice(y=(T(3.5), T(6.5)))
  d = convert(RectilinearGrid, cartgrid(10, 10))
  r, c = TB.apply(f, d)
  @test r isa RectilinearGrid
  @test r == convert(RectilinearGrid, CartesianGrid((10, 4), cart(0, 3), T.((1, 1))))
end

@testitem "Repair(0)" setup = [Setup] begin
  @test !isaffine(Repair)
  poly = PolyArea(cart.([(0, 0), (1, 0), (1, 0), (1, 1), (0, 1), (0, 1)]))
  rpoly = poly |> Repair(0)
  @test nvertices(rpoly) == 4
  @test vertices(rpoly) == cart.([(0, 0), (1, 0), (1, 1), (0, 1)])

  points = cart.([(0, 0), (1, 0), (0, 1), (1, 1), (1, 0), (2, 0), (1, 1), (2, 1)])
  connec = connect.([(1, 2, 4, 3), (5, 6, 8, 7)], Quadrangle)
  mesh = SimpleMesh(points, connec)
  rmesh = mesh |> Repair(0)
  @test nvertices(rmesh) == 6
  @test nelements(rmesh) == 2
  @test vertices(rmesh) == cart.([(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (2, 1)])
  @test collect(topology(rmesh)) == connect.([(1, 2, 4, 3), (2, 5, 6, 4)], Quadrangle)
end

@testitem "Repair(1)" setup = [Setup] begin
  # a tetrahedron with an unused vertex
  points = cart.([(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)])
  connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)])
  mesh = SimpleMesh(points, connec)
  rmesh = mesh |> Repair(1)
  @test nvertices(rmesh) == nvertices(mesh) - 1
  @test nelements(rmesh) == nelements(mesh)
  @test cart(5, 5, 5) ∉ vertices(rmesh)
end

@testitem "Repair(2)" setup = [Setup] begin end

@testitem "Repair(3)" setup = [Setup] begin end

@testitem "Repair(4)" setup = [Setup] begin end

@testitem "Repair(5)" setup = [Setup] begin end

@testitem "Repair(6)" setup = [Setup] begin end

@testitem "Repair(7)" setup = [Setup] begin
  # mesh with inconsistent orientation
  points = randpoint3(6)
  connec = connect.([(1, 2, 3), (3, 4, 2), (4, 3, 5), (6, 3, 1)])
  mesh = SimpleMesh(points, connec)
  rmesh = mesh |> Repair(7)
  topo = topology(mesh)
  rtopo = topology(rmesh)
  e = collect(elements(topo))
  n = collect(elements(rtopo))
  @test n[1] == e[1]
  @test n[2] != e[2]
  @test n[3] != e[3]
  @test n[4] != e[4]
end

@testitem "Repair(8)" setup = [Setup] begin
  poly =
    PolyArea(cart.([(0.0, 0.0), (0.5, -0.5), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 1.5), (0.0, 1.0), (-0.5, 0.5)]))
  rpoly = poly |> Repair(8)
  @test nvertices(rpoly) == 4
  @test vertices(rpoly) == cart.([(0.5, -0.5), (1.5, 0.5), (0.5, 1.5), (-0.5, 0.5)])

  # degenerate triangle with repeated vertices
  poly = PolyArea(cart.([(0, 0), (1, 1), (1, 1)]))
  rpoly = poly |> Repair(8)
  @test !hasholes(rpoly)
  @test rings(rpoly) == [Ring(cart(0, 0))]
  @test vertices(rpoly) == [cart(0, 0)]
end

@testitem "Repair(9)" setup = [Setup] begin
  quad = Quadrangle(cart(0, 1, 0), cart(1, 1, 0), cart(1, 0, 0), cart(0, 0, 0))
  repair = Repair(9)
  rquad, cache = TB.apply(repair, quad)
  @test rquad isa Quadrangle
  @test rquad == quad

  outer = Ring(cart(6, 4), cart(6, 7), cart(1, 6), cart(1, 1), cart(5, 2))
  inner1 = Ring(cart(3, 3), cart(3, 4), cart(4, 3))
  inner2 = Ring(cart(2, 5), cart(2, 6), cart(3, 5))
  poly = PolyArea([outer, inner1, inner2])
  repair = Repair(9)
  rpoly, cache = TB.apply(repair, poly)
  @test rpoly == PolyArea([outer, inner2, inner1])
end

@testitem "Repair(10)" setup = [Setup] begin
  outer = Ring(cart.([(0, 0), (0, 3), (2, 3), (2, 2), (3, 2), (3, 0)]))
  inner = Ring(cart.([(1, 1), (1, 2), (2, 2), (2, 1)]))
  poly = PolyArea(outer, inner)
  repair = Repair(10)
  rpoly, cache = TB.apply(repair, poly)
  @test nvertices(rpoly) == nvertices(poly)
  @test length(first(rings(rpoly))) > length(first(rings(poly)))
  opoly = TB.revert(repair, rpoly, cache)
  @test opoly ≈ poly
end

@testitem "Repair(11)" setup = [Setup] begin
  outer = cart.([(0, 0), (0, 2), (2, 2), (2, 0)])
  inner = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
  poly = PolyArea(outer, inner)
  repair = Repair(11)
  rpoly, cache = TB.apply(repair, poly)
  router, rinner = rings(rpoly)
  @test router == Ring(cart.([(0, 0), (2, 0), (2, 2), (0, 2)]))
  @test rinner == Ring(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))
end

@testitem "Repair(12)" setup = [Setup] begin
  poly = PolyArea(cart.([(0, 0), (1, 0)]))
  repair = Repair(12)
  rpoly, cache = TB.apply(repair, poly)
  @test rpoly == PolyArea(cart.([(0, 0), (0.5, 0.0), (1, 0)]))

  outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
  inner = cart.([(1, 2), (2, 3)])
  poly = PolyArea(outer, inner)
  repair = Repair(12)
  rpoly, cache = TB.apply(repair, poly)
  @test rpoly == PolyArea(outer)
end

@testitem "Repair fallbacks" setup = [Setup] begin
  quad = Quadrangle(cart(0, 1, 0), cart(1, 1, 0), cart(1, 0, 0), cart(0, 0, 0))
  repair = Repair(10)
  rquad, cache = TB.apply(repair, quad)
  @test rquad isa Quadrangle
  @test rquad == quad

  poly1 = PolyArea(cart.([(0, 0), (0, 2), (2, 2), (2, 0)]))
  poly2 = PolyArea(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))
  multi = Multi([poly1, poly2])
  repair = Repair(11)
  rmulti, cache = TB.apply(repair, multi)
  @test rmulti == Multi([repair(poly1), repair(poly2)])

  poly1 = PolyArea(cart.([(0, 0), (0, 2), (2, 2), (2, 0)]))
  poly2 = PolyArea(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))
  gset = GeometrySet([poly1, poly2])
  repair = Repair(11)
  rgset, cache = TB.apply(repair, gset)
  @test rgset == GeometrySet([repair(poly1), repair(poly2)])
end

@testitem "Repair IO" setup = [Setup] begin
  for K in 0:12
    repair = Repair(K)
    @test sprint(show, repair) == "Repair(K: $K)"
    @test sprint(show, MIME"text/plain"(), repair) == """
    Repair transform
    └─ K: $K"""
  end
end

@testitem "Bridge" setup = [Setup] begin
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

  # make sure that result is inferred
  @inferred poly |> Bridge(T(0.1))

  # unique and bridges
  poly = PolyArea(cart.([(0, 0), (1, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1), (0, 1)]))
  cpoly = poly |> Repair(0) |> Bridge()
  @test cpoly == PolyArea(cart.([(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)]))

  # basic ngon tests
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  @test (t |> Bridge() |> boundary) == boundary(t)
  q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  @test (q |> Bridge() |> boundary) == boundary(q)

  # bridges between holes
  outer = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  hole1 = Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))
  hole2 = Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))
  poly = PolyArea([outer, reverse(hole1), reverse(hole2)])
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

@testitem "Smoothing" setup = [Setup] begin
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
