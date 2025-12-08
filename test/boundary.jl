@testitem "boundary" setup = [Setup] begin
  # point
  @test isnothing(boundary(cart(1)))
  @test isnothing(boundary(cart(1, 2)))
  @test isnothing(boundary(cart(1, 2, 3)))
  @test embedboundary(cart(1)) == cart(1)
  @test embedboundary(cart(1, 2)) == cart(1, 2)
  @test embedboundary(cart(1, 2, 3)) == cart(1, 2, 3)

  # ray
  r = Ray(cart(0, 0), vector(1, 1))
  @test boundary(r) == cart(0, 0)
  @test embedboundary(r) == r

  # line
  l = Line(cart(0, 0), cart(1, 1))
  @test isnothing(boundary(l))
  @test embedboundary(l) == l

  # Bezier curve
  b = BezierCurve(cart(0, 0), cart(0.5, 1), cart(1, 0))
  @test boundary(b) == Multi([cart(0, 0), cart(1, 0)])
  @test embedboundary(b) == b
  b = BezierCurve(cart(0, 0), cart(1, 1))
  @test boundary(b) == Multi([cart(0, 0), cart(1, 1)])
  @test embedboundary(b) == b

  # parametrized curve
  c = ParametrizedCurve(t -> Point(Polar(T(1), T(t))), (T(0), T(2π)))
  @test isnothing(boundary(c))
  @test embedboundary(c) == c
  c = ParametrizedCurve(t -> cart(cospi(t), sinpi(t)), (T(0), T(1)))
  @test boundary(c) == Multi([cart(1, 0), cart(-1, 0)])

  # plane
  p = Plane(cart(0, 0, 0), vector(1, 0, 0), vector(0, 1, 0))
  @test isnothing(boundary(p))
  @test embedboundary(p) == p

  # box
  b = Box(cart(0), cart(1))
  @test boundary(b) == embedboundary(b) == Multi([cart(0), cart(1)])
  b = Box(cart(1, 2), cart(3, 4))
  v = cart.([(1, 2), (3, 2), (3, 4), (1, 4)])
  @test boundary(b) == embedboundary(b) == Ring(v)
  b = Box(cart(1, 2, 3), cart(4, 5, 6))
  v = cart.([(1, 2, 3), (4, 2, 3), (4, 5, 3), (1, 5, 3), (1, 2, 6), (4, 2, 6), (4, 5, 6), (1, 5, 6)])
  c = connect.([(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)])
  @test boundary(b) == embedboundary(b) == SimpleMesh(v, c)
  b = Box(cart(0, 0), cart(1, 1))
  @test boundary(b) == embedboundary(b) == Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  b = Box(latlon(0, 0), latlon(1, 1))
  @test boundary(b) == embedboundary(b) == Ring(latlon.([(0, 0), (0, 1), (1, 1), (1, 0)]))
  b = Box(cart(0, 0, 0), cart(1, 1, 1))
  m = boundary(b)
  @test m isa Mesh
  @test nvertices(m) == 8
  @test nelements(m) == 6
  b = Box(cart(0, 0, 0), cart(1, 1, 1))
  m = embedboundary(b)
  @test m isa Mesh
  @test nvertices(m) == 8
  @test nelements(m) == 6

  # disk
  p = Plane(cart(0, 0, 0), vector(0, 0, 1))
  d = Disk(p, T(2))
  @test boundary(d) == Circle(p, T(2))
  @test embedboundary(d) == d

  # circle
  p = Plane(cart(0, 0, 0), vector(0, 0, 1))
  c = Circle(p, T(2))
  @test isnothing(boundary(c))
  @test embedboundary(c) == c

  # ball
  b = Ball(cart(0, 0), T(1))
  @test boundary(b) == embedboundary(b) == Sphere(cart(0, 0), T(1))
  b = Ball(cart(0, 0, 0), T(1))
  @test boundary(b) == embedboundary(b) == Sphere(cart(0, 0, 0), T(1))

  # sphere
  s = Sphere(cart(0, 0), T(1))
  @test isnothing(boundary(s))
  @test embedboundary(s) == s
  s = Sphere(cart(0, 0, 0), T(1))
  @test isnothing(boundary(s))
  @test embedboundary(s) == s

  # ellipsoid
  e = Ellipsoid((T(3), T(2), T(1)))
  @test isnothing(boundary(e))
  @test embedboundary(e) == e

  # torus
  t = Torus(T.((1, 1, 1)), T.((1, 0, 0)), 2, 1)
  @test isnothing(boundary(t))
  @test embedboundary(t) == t

  # paraboloid surface
  p = ParaboloidSurface(cart(0, 0, 0), T(1), T(2))
  @test isnothing(boundary(p))
  @test embedboundary(p) == p

  # cylinder
  c = Cylinder(cart(0, 0, 0), cart(0, 0, 1), T(1))
  @test boundary(c) == embedboundary(c) == CylinderSurface(cart(0, 0, 0), cart(0, 0, 1), T(1))

  # cylinder surface
  c = CylinderSurface(T(2))
  @test isnothing(boundary(c))
  @test embedboundary(c) == c

  # cone
  p = Plane(cart(0, 0, 0), vector(0, 0, 1))
  d = Disk(p, T(2))
  a = cart(0, 0, 1)
  c = Cone(d, a)
  @test boundary(c) == embedboundary(c) == ConeSurface(d, a)

  # cone surface
  p = Plane(cart(0, 0, 0), vector(0, 0, 1))
  d = Disk(p, T(2))
  a = cart(0, 0, 1)
  s = ConeSurface(d, a)
  @test isnothing(boundary(s))
  @test embedboundary(s) == s

  # frustum
  pb = Plane(cart(0, 0, 0), vector(0, 0, 1))
  db = Disk(pb, T(1))
  pt = Plane(cart(0, 0, 10), vector(0, 0, 1))
  dt = Disk(pt, T(2))
  f = Frustum(db, dt)
  @test boundary(f) == embedboundary(f) == FrustumSurface(db, dt)

  # frustum surface
  pb = Plane(cart(0, 0, 0), vector(0, 0, 1))
  db = Disk(pb, T(1))
  pt = Plane(cart(0, 0, 10), vector(0, 0, 1))
  dt = Disk(pt, T(2))
  f = FrustumSurface(db, dt)
  @test isnothing(boundary(f))
  @test embedboundary(f) == f

  # segment
  s = Segment(cart(0), cart(1))
  @test boundary(s) == Multi([cart(0), cart(1)])
  @test embedboundary(s) == boundary(s)
  s = Segment(cart(0, 0), cart(1, 1))
  @test boundary(s) == Multi([cart(0, 0), cart(1, 1)])
  @test embedboundary(s) == s
  s = Segment(cart(0, 0, 0), cart(1, 1, 1))
  @test boundary(s) == Multi([cart(0, 0, 0), cart(1, 1, 1)])
  @test embedboundary(s) == s

  # chain
  c = Rope(cart.([(0,), (1,), (1,), (0,)]))
  @test boundary(c) == Multi(cart.([(0,), (0,)]))
  @test embedboundary(c) == boundary(c)
  c = Rope(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test boundary(c) == Multi(cart.([(0, 0), (0, 1)]))
  @test embedboundary(c) == c
  c = Rope(cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
  @test boundary(c) == Multi(cart.([(0, 0, 0), (0, 1, 0)]))
  @test embedboundary(c) == c
  c = Ring(cart.([(0,), (1,), (1,), (0,)]))
  @test isnothing(boundary(c))
  @test embedboundary(c) == boundary(c)
  c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test isnothing(boundary(c))
  @test embedboundary(c) == c
  c = Ring(cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
  @test isnothing(boundary(c))
  @test embedboundary(c) == c

  # triangle
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  @test boundary(t) == embedboundary(t) == first(rings(t))
  t = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0))
  @test boundary(t) == first(rings(t))
  @test embedboundary(t) == t

  # quadrangle
  q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  @test boundary(q) == embedboundary(q) == first(rings(q))
  q = Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))
  @test boundary(q) == first(rings(q))
  @test embedboundary(q) == q

  # polyarea
  outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
  hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
  hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
  poly = PolyArea([outer, hole1, hole2])
  @test boundary(poly) == Multi(rings(poly))
  @test embedboundary(poly) == boundary(poly)
  outer = cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
  hole1 = cart.([(0.2, 0.2, 0.0), (0.4, 0.2, 0.0), (0.4, 0.4, 0.0), (0.2, 0.4, 0.0)])
  hole2 = cart.([(0.6, 0.2, 0.0), (0.8, 0.2, 0.0), (0.8, 0.4, 0.0), (0.6, 0.4, 0.0)])
  poly = PolyArea([outer, hole1, hole2])
  @test boundary(poly) == Multi(rings(poly))
  @test embedboundary(poly) == poly

  # tetrahedron
  t = Tetrahedron(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1))
  @test boundary(t) == embedboundary(t)
  m = boundary(t)
  n = normal.(m)
  @test m isa Mesh
  @test nvertices(m) == 4
  @test nelements(m) == 4
  @test n[1] == vector(0, 0, -1)
  @test n[2] == vector(0, -1, 0)
  @test n[3] == vector(-1, 0, 0)
  @test all(>(T(0) * u"m"), n[4])

  # hexahedron
  h = Hexahedron(
    cart(0, 0, 0),
    cart(1, 0, 0),
    cart(1, 1, 0),
    cart(0, 1, 0),
    cart(0, 0, 1),
    cart(1, 0, 1),
    cart(1, 1, 1),
    cart(0, 1, 1)
  )
  @test boundary(h) == embedboundary(h)
  m = boundary(h)
  @test m isa Mesh
  @test nvertices(m) == 8
  @test nelements(m) == 6

  # pyramid
  p = Pyramid(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0), cart(0, 0, 1))
  @test boundary(p) == embedboundary(p)
  m = boundary(p)
  @test m isa Mesh
  @test nelements(m) == 5
  @test m[1] isa Quadrangle
  @test m[2] isa Triangle
  @test m[3] isa Triangle
  @test m[4] isa Triangle
  @test m[5] isa Triangle

  # wedge
  w = Wedge(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1), cart(1, 0, 1), cart(0, 1, 1))
  @test boundary(w) == embedboundary(w)
  m = boundary(w)
  @test m isa Mesh
  @test nelements(m) == 5
  @test m[1] isa Triangle
  @test m[2] isa Triangle
  @test m[3] isa Quadrangle
  @test m[4] isa Quadrangle
  @test m[5] isa Quadrangle

  # multi-geometry
  outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
  hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
  hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
  poly = PolyArea([outer, hole1, hole2])
  multi = Multi([poly, poly])
  @test boundary(multi) == merge(boundary(poly), boundary(poly))
  @test embedboundary(multi) == merge(embedboundary(poly), embedboundary(poly))
  box1 = Box(cart(0, 0), cart(1, 1))
  box2 = Box(cart(1, 1), cart(2, 2))
  mbox = Multi([box1, box2])
  mchn = boundary(mbox)
  noth = boundary(mchn)
  @test mchn isa Multi
  @test isnothing(noth)
  @test length(mchn) == T(8) * u"m"
  @test embedboundary(mbox) == boundary(mbox)
  @test embedboundary(mchn) == mchn
  box1 = Box(cart(0, 0, 0), cart(1, 1, 1))
  box2 = Box(cart(1, 1, 1), cart(2, 2, 2))
  mbox = Multi([box1, box2])
  mesh = boundary(mbox)
  @test mesh isa Mesh
  @test nvertices(mesh) == 16
  @test nelements(mesh) == 12
  @test embedboundary(mbox) == boundary(mbox)

  # transformed geometry
  b = Box(cart(0, 0), cart(1, 1))
  t = Translate(T(1), T(2))
  tb = TransformedGeometry(b, t)
  @test boundary(tb) == t(boundary(b))
  @test embedboundary(tb) == boundary(tb)
  b = Box(latlon(0, 0), latlon(1, 1))
  t = Proj(Mercator)
  tb = TransformedGeometry(b, t)
  @test boundary(tb) == t(boundary(b))
  @test embedboundary(tb) == boundary(tb)
  b = Box(cart(0, 0, 0), cart(1, 1, 1))
  t = Translate(T(1), T(2), T(3))
  tb = TransformedGeometry(b, t)
  @test boundary(tb) isa Mesh
  r = Ring(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  t = Translate(T(1), T(2))
  tr = TransformedGeometry(r, t)
  @test isnothing(boundary(tr))
  r = Ring(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  t = Translate(T(1), T(2))
  tr = TransformedGeometry(r, t)
  @test embedboundary(tr) == tr

  # transformed geometry with paramdim == 2 && embeddim == 3
  b = Box(cart(0, 0), cart(1, 1))
  t = ReinterpretCoords(Cartesian, LatLon)
  tb = TransformedGeometry(b, t)
  @test paramdim(tb) == 2
  @test embeddim(tb) == 3
  @test embedboundary(tb) == tb
end

@testitem "boundarypoints" setup = [Setup] begin
  p = cart(0, 0)
  @test boundarypoints(p) == [cart(0, 0)]

  sphere = Sphere(cart(0, 0), T(1))
  points = boundarypoints(sphere)
  @test all(∈(sphere), points)

  ball = Ball(cart(0, 0), T(1))
  points = boundarypoints(ball)
  @test all(∈(boundary(ball)), points)

  verts = [cart(0, 0), cart(1, 1)]
  segment = Segment(verts...)
  points = boundarypoints(segment)
  @test points == verts

  verts = [cart(0, 0), cart(1, 0), cart(1, 1)]
  triangle = Triangle(verts...)
  points = boundarypoints(triangle)
  @test points == verts

  verts = [cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1)]
  quadrangle = Quadrangle(verts...)
  points = boundarypoints(quadrangle)
  @test points == verts

  tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  multi = Multi([tri, quad])
  points = boundarypoints(multi)
  @test points == [boundarypoints(tri); boundarypoints(quad)]

  box = Box(cart(0, 0), cart(1, 1))
  trans = Translate(T(1), T(2))
  tbox = TransformedGeometry(box, trans)
  points = boundarypoints(tbox)
  @test points == trans.(boundarypoints(box))

  box = Box(latlon(0, 0), latlon(45, 45))
  trans = Proj(Mercator)
  tbox = TransformedGeometry(box, trans)
  points = boundarypoints(tbox)
  @test points == vertices(discretize(boundary(tbox)))
end
