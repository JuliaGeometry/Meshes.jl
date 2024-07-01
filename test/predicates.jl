@testset "Predicates" begin
  @testset "issimplex" begin
    @test issimplex(Segment)
    @test issimplex(Segment(cart(0, 0), cart(1, 0)))

    @test issimplex(Triangle)
    @test issimplex(Triangle(cart(0, 0), cart(1, 0), cart(0, 1)))

    @test issimplex(Tetrahedron)
    @test issimplex(Tetrahedron(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1)))
  end

  @testset "isconvex" begin
    # primitives
    r = Ray(cart(0, 0), vector(1, 1))
    @test isconvex(r)
    l = Line(cart(0, 0), cart(1, 1))
    @test isconvex(l)
    p = Plane(cart(0, 0, 0), vector(1, 0, 0), vector(0, 1, 0))
    @test isconvex(p)
    b = Box(cart(0), cart(1))
    @test isconvex(b)
    b = Box(cart(0, 0), cart(1, 1))
    @test isconvex(b)
    b = Box(cart(0, 0, 0), cart(1, 1, 1))
    b = Ball(cart(1, 2, 3), T(5))
    @test isconvex(b)
    @test isconvex(b)
    s = Sphere(cart(0, 0), T(1))
    @test !isconvex(s)
    s = Sphere(cart(0, 0, 0), T(1))
    @test !isconvex(s)
    d = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    @test isconvex(d)
    c = Circle(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    @test !isconvex(c)
    b = BezierCurve(cart.([(0, 0), (1, 0), (2, 0)]))
    @test isconvex(b)
    b = BezierCurve(cart.([(0, 0), (1, 1), (2, 2)]))
    @test isconvex(b)
    b = BezierCurve(cart.([(0, 0)]))
    @test isconvex(b)
    b = BezierCurve(cart.([(0, 0), (1, 0)]))
    @test isconvex(b)
    b = BezierCurve(cart.([(0, 0), (5, 3), (-10, 3), (17, 20)]))
    @test !isconvex(b)
    b = BezierCurve(cart.([(5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (5, 11)]))
    @test isconvex(b)
    P = typeof(cart(0, 0))
    b = BezierCurve(P[])
    @test isconvex(b)
    c = Cylinder(Plane(cart(1, 2, 3), vector(0, 0, 1)), Plane(cart(4, 5, 6), vector(0, 0, 1)), T(5))
    @test isconvex(c)
    c = CylinderSurface(T(2))
    @test !isconvex(c)
    d = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    a = cart(0, 0, 1)
    c = Cone(d, a)
    @test isconvex(c)
    d = Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(2))
    a = cart(0, 0, 1)
    c = ConeSurface(d, a)
    @test !isconvex(c)
    t = Torus(T.((1, 1, 1)), T.((1, 0, 0)), 2, 1)
    @test !isconvex(t)

    # polytopes
    s = Segment(cart(0, 0), cart(1, 1))
    @test isconvex(s)
    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test isconvex(t)
    q1 = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    q2 = Quadrangle(cart(0.8, 0.8), cart(1, 0), cart(1, 1), cart(0, 1))
    q3 = Quadrangle(cart(0, 0), cart(0.2, 0.8), cart(1, 1), cart(0, 1))
    q4 = Quadrangle(cart(0, 0), cart(1, 0), cart(0.2, 0.2), cart(0, 1))
    q5 = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0.8, 0.2))
    @test isconvex(q1)
    @test !isconvex(q2)
    @test !isconvex(q3)
    @test !isconvex(q4)
    @test !isconvex(q5)
    q1 = Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))
    q2 = Quadrangle(cart(0.8, 0.8, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))
    q3 = Quadrangle(cart(0, 0, 0), cart(0.2, 0.8, 0), cart(1, 1, 0), cart(0, 1, 0))
    q4 = Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(0.2, 0.2, 0), cart(0, 1, 0))
    q5 = Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0.8, 0.2, 0))
    @test isconvex(q1)
    @test !isconvex(q2)
    @test !isconvex(q3)
    @test !isconvex(q4)
    @test !isconvex(q5)
    t = Tetrahedron(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1))
    @test isconvex(t)
    outer = cart.([(6, 1), (2, 10), (10, 16), (18, 10), (14, 1)])
    inner = cart.([(5, 7), (10, 12), (15, 7)])
    pent = Pentagon(outer...)
    tri = Triangle(inner...)
    poly = PolyArea([outer, inner])
    multi = Multi([poly, tri])
    @test isconvex(pent)
    @test isconvex(tri)
    @test !isconvex(poly)
    @test isconvex(multi)
    outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly1 = PolyArea(outer)
    poly2 = PolyArea([outer, hole1, hole2])
    @test isconvex(poly1)
    @test !isconvex(poly2)
    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0.5, 0.5), (0, 1)]))
    @test !isconvex(poly)
  end

  @testset "isparametrized" begin
    # primitives
    @test isparametrized(Ray)
    @test isparametrized(Line)
    @test isparametrized(Plane)
    @test isparametrized(Box)
    @test isparametrized(Ball)
    @test isparametrized(Sphere)
    @test isparametrized(Ellipsoid)
    @test isparametrized(Disk)
    @test isparametrized(Circle)
    @test isparametrized(BezierCurve)
    @test isparametrized(Cylinder)
    @test isparametrized(CylinderSurface)
    @test isparametrized(ConeSurface)
    @test isparametrized(ParaboloidSurface)
    @test isparametrized(Torus)

    # polytopes
    @test isparametrized(Segment)
    @test isparametrized(Triangle)
    @test isparametrized(Quadrangle)
    @test isparametrized(Hexahedron)
  end

  @testset "isperiodic" begin
    # primitives
    @test isperiodic(Box{1}) == (false,)
    @test isperiodic(Box{2}) == (false, false)
    @test isperiodic(Box{3}) == (false, false, false)
    @test isperiodic(Ball{2}) == (false, true)
    @test isperiodic(Ball{3}) == (false, true, true)
    @test isperiodic(Sphere{2}) == (true,)
    @test isperiodic(Sphere{3}) == (true, true)
    @test isperiodic(Ellipsoid) == (true, true)
    @test isperiodic(ParaboloidSurface) == (false, true)
    @test isperiodic(Torus) == (true, true)

    # polytopes
    @test isperiodic(Segment) == (false,)
    @test isperiodic(Quadrangle) == (false, false)
    @test isperiodic(Hexahedron) == (false, false, false)

    @test isperiodic(cartgrid(10, 10)) == (false, false)
    @test isperiodic(cartgrid(10, 10, 10)) == (false, false, false)
  end

  @testset "in" begin
    h = first(cartgrid(10, 10, 10))
    @test cart(0, 0, 0) ∈ h
    @test cart(0.5, 0.5, 0.5) ∈ h
    @test cart(-1, 0, 0) ∉ h
    @test cart(0, 2, 0) ∉ h

    outer = [latlon(0, 0), latlon(1, 0), latlon(1, 1), latlon(0, 1)]
    hole1 = [latlon(0.2, 0.2), latlon(0.4, 0.2), latlon(0.4, 0.4), latlon(0.2, 0.4)]
    hole2 = [latlon(0.6, 0.2), latlon(0.8, 0.2), latlon(0.8, 0.4), latlon(0.6, 0.4)]
    poly = PolyArea([outer, hole1, hole2])
    @test all(p ∈ poly for p in outer)
    @test latlon(0.5, 0.5) ∈ poly
    @test latlon(0.2, 0.6) ∈ poly
    @test latlon(1.5, 0.5) ∉ poly
    @test latlon(-0.5, 0.5) ∉ poly
    @test latlon(0.25, 0.25) ∉ poly
    @test latlon(0.75, 0.25) ∉ poly
    @test latlon(0.75, 0.75) ∈ poly
  end

  @testset "issubset" begin
    p = cart(0.5, 0.5)
    box = Box(cart(0, 0), cart(1, 1))
    ball = Ball(cart(0, 0))
    tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test p ⊆ box
    @test p ⊆ ball
    @test p ⊆ tri
    @test p ⊆ quad
    @test p ⊆ p
    @test quad ⊆ quad

    s1 = Segment(cart(0, 0), cart(1, 1))
    s2 = Segment(cart(0.5, 0.5), cart(1, 1))
    s3 = Segment(cart(0, 0), cart(0.5, 0.5))
    @test s2 ⊆ s1
    @test s3 ⊆ s1
    @test s1 ⊆ s1

    seg = Segment(cart(0, 0), cart(1, 1))
    box = Box(cart(0, 0), cart(1, 1))
    ball = Ball(cart(0, 0))
    @test seg ⊆ box
    @test !(seg ⊆ ball)

    t1 = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    t2 = Triangle(cart(0, 0), cart(1, 0), cart(0.8, 0.8))
    t3 = Triangle(cart(0, 0), cart(1, 0), cart(1.1, 1.1))
    @test t2 ⊆ t1
    @test !(t3 ⊆ t1)

    tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    box = Box(cart(0, 0), cart(1, 1))
    ball = Ball(cart(0, 0))
    quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    pent = Pentagon(cart(0, 0), cart(1, 0), cart(1, 1), cart(0.5, 1.5), cart(0, 1))
    poly = PolyArea(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test tri ⊆ quad
    @test !(quad ⊆ tri)
    @test tri ⊆ box
    @test !(box ⊆ tri)
    @test !(tri ⊆ ball)
    @test !(ball ⊆ tri)
    @test tri ⊆ pent
    @test !(pent ⊆ tri)
    @test quad ⊆ pent
    @test !(pent ⊆ quad)
    @test tri ⊆ poly
    @test !(poly ⊆ tri)
    @test quad ⊆ poly
    @test poly ⊆ quad

    quad1 = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    quad2 = Quadrangle(cart(0, 0), cart(1.1, 0), cart(1, 1), cart(0, 1))
    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    multi = Multi([poly])
    @test quad1 ⊆ poly
    @test !(quad2 ⊆ poly)
    @test quad1 ⊆ multi
    @test !(quad2 ⊆ multi)

    p1 = cart(-1.0, 0.0)
    p2 = cart(0.0, 0.0)
    p3 = cart(1.0, 0.0)
    l1 = Line(p1, p3)
    l2 = Line(p2, p3)
    @test l1 ⊆ l2
    @test l2 ⊆ l1
    @test l1 ⊆ l1
    @test l2 ⊆ l2

    pts1 = cart.([(5, 7), (10, 12), (15, 7)])
    pts2 = cart.([(6, 1), (2, 10), (10, 16), (18, 10), (14, 1)])
    pent = Pentagon(pts2...)
    tri = Triangle(pts1...)
    poly1 = PolyArea(pts2)
    poly2 = PolyArea([pts2, pts1])
    multi = Multi([poly2, tri])
    @test tri ⊆ pent
    @test tri ⊆ poly1
    @test tri ⊈ poly2
    @test tri ⊆ multi
    @test pent ⊆ poly1
    @test pent ⊈ poly2
    @test pent ⊆ multi

    poly1 = PolyArea(cart.([(4, 12), (11, 11), (16, 8), (16, 1), (13, -2), (2, -2), (-3, 4), (-2, 8)]))
    poly2 = PolyArea(cart.([(3, 0), (1, 2), (3, 4), (1, 6), (4, 7), (10, 7), (11, 4), (9, 0)]))
    poly3 = PolyArea(cart.([(3, 2), (4, 4), (3, 8), (12, 8), (14, 4), (12, 1)]))
    poly4 = PolyArea(cart.([(8, 2), (5, 4), (5, 6), (9, 6), (10, 4)]))
    poly5 = PolyArea(cart.([(3, 9), (6, 11), (10, 10), (10, 9)]))
    @test poly2 ⊆ poly1
    @test poly3 ⊆ poly1
    @test poly4 ⊆ poly1
    @test poly5 ⊆ poly1
    @test poly4 ⊆ poly2
    @test poly4 ⊆ poly3
    @test poly5 ⊈ poly2
    @test poly5 ⊈ poly3
  end

  @testset "intersects" begin
    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    q = Quadrangle(cart(1, 1), cart(2, 1), cart(2, 2), cart(1, 2))
    @test intersects(t, t)
    @test intersects(q, q)
    @test !intersects(t, q)
    @test !intersects(q, t)

    t = Triangle(cart(1, 0), cart(2, 0), cart(1, 1))
    q = Quadrangle(cart(1.3, 0.5), cart(2.3, 0.5), cart(2.3, 1.5), cart(1.3, 1.5))
    @test intersects(t, t)
    @test intersects(q, q)
    @test intersects(t, q)
    @test intersects(q, t)

    t = Triangle(cart(1, 0), cart(2, 0), cart(1, 1))
    q = Quadrangle(cart(1.3, 0.5), cart(2.3, 0.5), cart(2.3, 1.5), cart(1.3, 1.5))
    m = Multi([t, q])
    @test intersects(m, t)
    @test intersects(t, m)
    @test intersects(m, q)
    @test intersects(q, m)
    @test intersects(m, m)

    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    b = Ball(cart(0, 0), T(1))
    @test intersects(t, t)
    @test intersects(b, b)
    @test intersects(t, b)
    @test intersects(b, t)

    t = Triangle(cart(1, 0), cart(2, 0), cart(1, 1))
    b = Ball(cart(0, 0), T(1))
    @test intersects(t, t)
    @test intersects(b, b)
    @test intersects(t, b)
    @test intersects(b, t)

    t = Triangle(cart(1, 0), cart(2, 0), cart(1, 1))
    b = Ball(cart(-0.01, 0), T(1))
    @test intersects(t, t)
    @test intersects(b, b)
    @test !intersects(t, b)
    @test !intersects(b, t)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/250
    t1 = Triangle(cart(0, 0, 0), cart(2, 0, 0), cart(1, 2, 0))
    t2 = Triangle(cart(1, 0, 0), cart(3, 0, 0), cart(2, 2, 0))
    t3 = Triangle(cart(3, 0, 0), cart(5, 0, 0), cart(4, 2, 0))
    @test intersects(t1, t2)
    @test intersects(t2, t3)
    @test !intersects(t1, t3)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/639
    r = Ray(cart(0.41169768366272996, 0.8990554132423699), vector(0.47249211625247445, 0.2523149692768657))
    b = Box(cart(1.0, 1.0), cart(5.0, 2.0))
    @test intersects(r, b)
    @test intersects(b, r)

    t = Triangle(cart(0, 0, 0), cart(2, 0, 0), cart(1, 2, 0))
    r1 = Ray(cart(1, 1, 1), vector(0, 0, -1))
    r2 = Ray(cart(1, 1, 1), vector(0, 0, 1))
    @test intersects(r1, t)
    @test intersects(t, r1)
    @test !intersects(r2, t)
    @test !intersects(t, r2)

    r = Ray(cart(0, 0), vector(1, 0))
    s1 = Sphere(cart(3, 0), T(1))
    s2 = Sphere(cart(0, 3), T(1))
    @test intersects(r, s1)
    @test !intersects(r, s2)

    # result doesn't change under translation
    t1 = Translate(T(10), T(0))
    t2 = Translate(T(0), T(10))
    t3 = Translate(T(-10), T(0))
    t4 = Translate(T(0), T(-10))
    for t in [t1, t2, t3, t4]
      @test intersects(t(r), t(s1))
      @test !intersects(t(r), t(s2))
    end

    # result doesn't change under rotation
    r1 = Rotate(Angle2d(T(π / 2)))
    r2 = Rotate(Angle2d(T(-π / 2)))
    r3 = Rotate(Angle2d(T(π)))
    r4 = Rotate(Angle2d(T(-π)))
    for t in [r1, r2, r3, r4]
      @test intersects(t(r), t(s1))
      @test !intersects(t(r), t(s2))
    end

    r = Ray(cart(0, 0), vector(1, 0))
    s = Sphere(cart(floatmax(Float32) / 2, 0), 1)
    @test intersects(r, s)

    r = Ray(cart(0, 0, 0), vector(1, 0, 0))
    s1 = Sphere(cart(5, 0, 1 - eps(T(1))), T(1))
    s2 = Sphere(cart(5, 0, 1 + eps(T(1))), T(1))
    @test intersects(r, s1)
    @test !intersects(r, s2)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/635
    q1 = Quadrangle(cart(4.0, 4.0, 0.0), cart(3.0, 3.0, 2.0), cart(3.0, 1.0, 2.0), cart(4.0, 0.0, 0.0))
    q2 = Quadrangle(cart(3.6, 3.0, 1.0), cart(5.6, 3.0, 1.0), cart(5.6, 1.0, 1.0), cart(3.6, 1.0, 1.0))
    q3 = Quadrangle(cart(3.6, 1.0, 1.0), cart(5.6, 1.0, 1.0), cart(5.6, -1.0, 1.0), cart(3.6, -1.0, 1.0))
    q4 = Quadrangle(cart(2.1, 1.0, 1.0), cart(4.1, 1.0, 1.0), cart(4.1, -1.0, 1.0), cart(2.1, -1.0, 1.0))
    @test !intersects(q1, q2)
    @test !intersects(q1, q3)
    @test intersects(q1, q1)
    @test intersects(q1, q4)

    h1 = Tetrahedron(cart(1, 1, 0), cart(4, 4, 0), cart(2.5, 2.5, 1.5), cart(1, 3, 2))
    h2 = Tetrahedron(cart(-1.0, 2.0, 1.0), cart(2.0, 1.0, 1.0), cart(-1.0, 4.0, 0.0), cart(0.5, 2.5, 1.5))
    h3 = Tetrahedron(cart(-1.3, 2.0, 1.0), cart(1.7, 1.0, 1.0), cart(-1.3, 4.0, 0.0), cart(0.2, 2.5, 1.5))
    @test intersects(h1, h2)
    @test !intersects(h1, h3)

    outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly1 = PolyArea(outer)
    poly2 = PolyArea([outer, hole1, hole2])
    ball1 = Ball(cart(0.5, 0.5), T(0.05))
    ball2 = Ball(cart(0.3, 0.3), T(0.05))
    ball3 = Ball(cart(0.7, 0.3), T(0.05))
    ball4 = Ball(cart(0.3, 0.3), T(0.15))
    @test intersects(poly1, poly1)
    @test intersects(poly2, poly2)
    @test intersects(poly1, poly2)
    @test intersects(poly2, poly1)
    @test intersects(poly1, ball1)
    @test intersects(poly2, ball1)
    @test intersects(poly1, ball2)
    @test !intersects(poly2, ball2)
    @test intersects(poly1, ball3)
    @test !intersects(poly2, ball3)
    @test intersects(poly1, ball4)
    @test intersects(poly2, ball4)
    mesh1 = discretize(poly1, DehnTriangulation())
    mesh2 = discretize(poly2, DehnTriangulation())
    @test intersects(mesh1, mesh1)
    @test intersects(mesh2, mesh2)
    @test intersects(mesh1, mesh2)
    @test intersects(mesh2, mesh1)

    p = cart(0.5, 0.5)
    ball = Ball(cart(0, 0), T(1))
    @test intersects(p, ball)
    @test intersects(ball, p)
    @test intersects(p, p)
    @test !intersects(p, p + vector(1, 1))

    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    box = Box(cart(0, 0), cart(2, 2))
    @test intersects(poly, box)

    b1 = Box(cart(0, 0), cart(2, 2))
    b2 = Box(cart(2, 0), cart(4, 2))
    p1 = cart(1, 1)
    p2 = cart(3, 1)
    m = Multi([b1, b2])
    @test intersects(p1, b1)
    @test !intersects(p2, b1)
    @test intersects(p2, b2)
    @test !intersects(p1, b2)
    @test intersects(m, p1)
    @test intersects(p1, m)
    @test intersects(m, p2)
    @test intersects(p2, m)

    s1 = Segment(cart(0, 0), cart(4, 4))
    s2 = Segment(cart(4, 0), cart(0, 4))
    s3 = Segment(cart(2, 0), cart(4, 2))
    @test intersects(s1, s2)
    @test intersects(s2, s3)
    @test !intersects(s1, s3)

    s1 = Segment(cart(4, 0), cart(0, 4))
    s2 = Segment(cart(4, 0), cart(8, 4))
    s3 = Segment(cart(0, 8), cart(8, 8))
    r1 = Rope(cart.([(0, 0), (4, 4), (8, 0)]))
    r2 = Ring(cart.([(0, 2), (4, 6), (8, 2)]))
    @test intersects(s1, r1)
    @test intersects(s2, r1)
    @test !intersects(s3, r1)
    @test intersects(s1, r2)
    @test intersects(s2, r2)
    @test !intersects(s3, r2)
    @test intersects(r1, r2)

    r1 = Rope(cart.([(0, 0), (2, 2), (4, 0)]))
    r2 = Rope(cart.([(3, 0), (5, 2), (7, 0)]))
    r3 = Rope(cart.([(6, 0), (8, 2), (10, 0)]))
    @test intersects(r1, r2)
    @test intersects(r2, r3)
    @test !intersects(r1, r3)

    r1 = Ring(cart.([(0, 0), (2, 2), (4, 0)]))
    r2 = Ring(cart.([(3, 0), (5, 2), (7, 0)]))
    r3 = Ring(cart.([(6, 0), (8, 2), (10, 0)]))
    @test intersects(r1, r2)
    @test intersects(r2, r3)
    @test !intersects(r1, r3)

    t = Triangle(cart(3, 1), cart(7, 5), cart(11, 1))
    q = Quadrangle(cart(2, 0), cart(2, 7), cart(12, 7), cart(12, 0))
    b = Box(cart(2, 0), cart(12, 7))
    s1 = Segment(cart(5, 2), cart(9, 2))
    s2 = Segment(cart(0, 3), cart(5, 3))
    s3 = Segment(cart(4, 4), cart(10, 4))
    s4 = Segment(cart(1, 6), cart(13, 6))
    s5 = Segment(cart(0, 9), cart(14, 9))
    r1 = Ring(cart.([(1, 2), (7, 8), (13, 2)]))
    r2 = Rope(cart.([(1, 2), (7, 8), (13, 2)]))
    @test intersects(s1, t)
    @test intersects(s2, t)
    @test intersects(s3, t)
    @test !intersects(s4, t)
    @test !intersects(s5, t)
    @test intersects(s1, q)
    @test intersects(s2, q)
    @test intersects(s3, q)
    @test intersects(s4, q)
    @test !intersects(s5, q)
    @test intersects(s1, b)
    @test intersects(s2, b)
    @test intersects(s3, b)
    @test intersects(s4, b)
    @test !intersects(s5, b)
    @test intersects(r1, t)
    @test !intersects(r2, t)
    @test intersects(r1, q)
    @test intersects(r2, q)
    @test intersects(r1, b)
    @test intersects(r2, b)

    # performance test
    b1 = Box(cart(0, 0), cart(3, 3))
    b2 = Box(cart(2, 2), cart(5, 5))
    @test intersects(b1, b2)
    @test intersects(b2, b1)
    @test @elapsed(intersects(b1, b2)) < 5e-5
    @test @allocated(intersects(b1, b2)) < 100

    # partial application
    points = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly = PolyArea(points)
    box = Box(cart(0, 0), cart(2, 2))
    @test intersects(box)(poly)
    @test all(intersects(box), points)

    # method ambiguities
    p = cart(3, 1)
    ring = Ring(cart.([(0, 0), (2, 2), (4, 0)]))
    rope = Rope(cart.([(2, 0), (4, 2), (6, 0)]))
    seg = Segment(cart(0, 1), cart(6, 1))
    multi = Multi([ring])
    @test intersects(p, ring)
    @test intersects(p, rope)
    @test intersects(p, seg)
    @test intersects(p, multi)
    @test intersects(ring, multi)
    @test intersects(rope, multi)
    @test intersects(seg, multi)
  end

  @testset "iscollinear" begin
    @test iscollinear(cart(0, 0), cart(1, 1), cart(2, 2))
  end

  @testset "iscoplanar" begin
    @test iscoplanar(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))
  end
end
