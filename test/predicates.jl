@testset "Predicates" begin
  @testset "issimplex" begin
    @test issimplex(Segment)
    @test issimplex(Segment(P2(0, 0), P2(1, 0)))

    @test issimplex(Triangle)
    @test issimplex(Triangle(P2(0, 0), P2(1, 0), P2(0, 1)))

    @test issimplex(Tetrahedron)
    @test issimplex(Tetrahedron(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0), P3(0, 0, 1)))
  end

  @testset "isconvex" begin
    # primitives
    r = Ray(P2(0, 0), V2(1, 1))
    @test isconvex(r)
    l = Line(P2(0, 0), P2(1, 1))
    @test isconvex(l)
    p = Plane(P3(0, 0, 0), V3(1, 0, 0), V3(0, 1, 0))
    @test isconvex(p)
    b = Box(P1(0), P1(1))
    @test isconvex(b)
    b = Box(P2(0, 0), P2(1, 1))
    @test isconvex(b)
    b = Box(P3(0, 0, 0), P3(1, 1, 1))
    b = Ball(P3(1, 2, 3), T(5))
    @test isconvex(b)
    @test isconvex(b)
    s = Sphere(P2(0, 0), T(1))
    @test !isconvex(s)
    s = Sphere(P3(0, 0, 0), T(1))
    @test !isconvex(s)
    d = Disk(Plane(P3(0, 0, 0), V3(0, 0, 1)), T(2))
    @test isconvex(d)
    c = Circle(Plane(P3(0, 0, 0), V3(0, 0, 1)), T(2))
    @test !isconvex(c)
    b = BezierCurve(P2[(0, 0), (1, 0), (2, 0)])
    @test isconvex(b)
    b = BezierCurve(P2[(0, 0), (1, 1), (2, 2)])
    @test isconvex(b)
    b = BezierCurve(P2[(0, 0)])
    @test isconvex(b)
    b = BezierCurve(P2[(0, 0), (1, 0)])
    @test isconvex(b)
    b = BezierCurve(P2[(0, 0), (5, 3), (-10, 3), (17, 20)])
    @test !isconvex(b)
    b = BezierCurve(P2[(5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (5, 11)])
    @test isconvex(b)
    b = BezierCurve(P2[])
    @test isconvex(b)
    c = Cylinder(Plane(P3(1, 2, 3), V3(0, 0, 1)), Plane(P3(4, 5, 6), V3(0, 0, 1)), T(5))
    @test isconvex(c)
    c = CylinderSurface(T(2))
    @test !isconvex(c)
    d = Disk(Plane(P3(0, 0, 0), V3(0, 0, 1)), T(2))
    a = P3(0, 0, 1)
    c = Cone(d, a)
    @test isconvex(c)
    d = Disk(Plane(P3(0, 0, 0), V3(0, 0, 1)), T(2))
    a = P3(0, 0, 1)
    c = ConeSurface(d, a)
    @test !isconvex(c)
    t = Torus(T.((1, 1, 1)), T.((1, 0, 0)), 2, 1)
    @test !isconvex(t)

    # polytopes
    s = Segment(P2(0, 0), P2(1, 1))
    @test isconvex(s)
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test isconvex(t)
    q1 = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    q2 = Quadrangle(P2(0.8, 0.8), P2(1, 0), P2(1, 1), P2(0, 1))
    q3 = Quadrangle(P2(0, 0), P2(0.2, 0.8), P2(1, 1), P2(0, 1))
    q4 = Quadrangle(P2(0, 0), P2(1, 0), P2(0.2, 0.2), P2(0, 1))
    q5 = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0.8, 0.2))
    @test isconvex(q1)
    @test !isconvex(q2)
    @test !isconvex(q3)
    @test !isconvex(q4)
    @test !isconvex(q5)
    q1 = Quadrangle(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0))
    q2 = Quadrangle(P3(0.8, 0.8, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0))
    q3 = Quadrangle(P3(0, 0, 0), P3(0.2, 0.8, 0), P3(1, 1, 0), P3(0, 1, 0))
    q4 = Quadrangle(P3(0, 0, 0), P3(1, 0, 0), P3(0.2, 0.2, 0), P3(0, 1, 0))
    q5 = Quadrangle(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0.8, 0.2, 0))
    @test isconvex(q1)
    @test !isconvex(q2)
    @test !isconvex(q3)
    @test !isconvex(q4)
    @test !isconvex(q5)
    t = Tetrahedron(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0), P3(0, 0, 1))
    @test isconvex(t)
    outer = P2[(6, 1), (2, 10), (10, 16), (18, 10), (14, 1)]
    inner = P2[(5, 7), (10, 12), (15, 7)]
    pent = Pentagon(outer...)
    tri = Triangle(inner...)
    poly = PolyArea([outer, inner])
    multi = Multi([poly, tri])
    @test isconvex(pent)
    @test isconvex(tri)
    @test !isconvex(poly)
    @test isconvex(multi)
    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly1 = PolyArea(outer)
    poly2 = PolyArea([outer, hole1, hole2])
    @test isconvex(poly1)
    @test !isconvex(poly2)
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0.5, 0.5), (0, 1)])
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

    @test isperiodic(CartesianGrid{T}(10, 10)) == (false, false)
    @test isperiodic(CartesianGrid{T}(10, 10, 10)) == (false, false, false)
  end

  @testset "in" begin
    h = first(CartesianGrid{T}(10, 10, 10))
    @test P3(0, 0, 0) ∈ h
    @test P3(0.5, 0.5, 0.5) ∈ h
    @test P3(-1, 0, 0) ∉ h
    @test P3(0, 2, 0) ∉ h
  end

  @testset "issubset" begin
    point = P2(0.5, 0.5)
    box = Box(P2(0, 0), P2(1, 1))
    ball = Ball(P2(0, 0))
    tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test point ⊆ box
    @test point ⊆ ball
    @test point ⊆ tri
    @test point ⊆ quad
    @test point ⊆ point
    @test quad ⊆ quad

    s1 = Segment(P2(0, 0), P2(1, 1))
    s2 = Segment(P2(0.5, 0.5), P2(1, 1))
    s3 = Segment(P2(0, 0), P2(0.5, 0.5))
    @test s2 ⊆ s1
    @test s3 ⊆ s1
    @test s1 ⊆ s1

    seg = Segment(P2(0, 0), P2(1, 1))
    box = Box(P2(0, 0), P2(1, 1))
    ball = Ball(P2(0, 0))
    @test seg ⊆ box
    @test !(seg ⊆ ball)

    t1 = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    t2 = Triangle(P2(0, 0), P2(1, 0), P2(0.8, 0.8))
    t3 = Triangle(P2(0, 0), P2(1, 0), P2(1.1, 1.1))
    @test t2 ⊆ t1
    @test !(t3 ⊆ t1)

    tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    box = Box(P2(0, 0), P2(1, 1))
    ball = Ball(P2(0, 0))
    quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    pent = Pentagon(P2(0, 0), P2(1, 0), P2(1, 1), P2(0.5, 1.5), P2(0, 1))
    poly = PolyArea(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
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

    quad1 = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    quad2 = Quadrangle(P2(0, 0), P2(1.1, 0), P2(1, 1), P2(0, 1))
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    multi = Multi([poly])
    @test quad1 ⊆ poly
    @test !(quad2 ⊆ poly)
    @test quad1 ⊆ multi
    @test !(quad2 ⊆ multi)

    p1 = P2(-1.0, 0.0)
    p2 = P2(0.0, 0.0)
    p3 = P2(1.0, 0.0)
    l1 = Line(p1, p3)
    l2 = Line(p2, p3)
    @test l1 ⊆ l2
    @test l2 ⊆ l1
    @test l1 ⊆ l1
    @test l2 ⊆ l2

    pts1 = P2[(5, 7), (10, 12), (15, 7)]
    pts2 = P2[(6, 1), (2, 10), (10, 16), (18, 10), (14, 1)]
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

    poly1 = PolyArea(P2[(4, 12), (11, 11), (16, 8), (16, 1), (13, -2), (2, -2), (-3, 4), (-2, 8)])
    poly2 = PolyArea(P2[(3, 0), (1, 2), (3, 4), (1, 6), (4, 7), (10, 7), (11, 4), (9, 0)])
    poly3 = PolyArea(P2[(3, 2), (4, 4), (3, 8), (12, 8), (14, 4), (12, 1)])
    poly4 = PolyArea(P2[(8, 2), (5, 4), (5, 6), (9, 6), (10, 4)])
    poly5 = PolyArea(P2[(3, 9), (6, 11), (10, 10), (10, 9)])
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
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    q = Quadrangle(P2(1, 1), P2(2, 1), P2(2, 2), P2(1, 2))
    @test intersects(t, t)
    @test intersects(q, q)
    @test !intersects(t, q)
    @test !intersects(q, t)

    t = Triangle(P2(1, 0), P2(2, 0), P2(1, 1))
    q = Quadrangle(P2(1.3, 0.5), P2(2.3, 0.5), P2(2.3, 1.5), P2(1.3, 1.5))
    @test intersects(t, t)
    @test intersects(q, q)
    @test intersects(t, q)
    @test intersects(q, t)

    t = Triangle(P2(1, 0), P2(2, 0), P2(1, 1))
    q = Quadrangle(P2(1.3, 0.5), P2(2.3, 0.5), P2(2.3, 1.5), P2(1.3, 1.5))
    m = Multi([t, q])
    @test intersects(m, t)
    @test intersects(t, m)
    @test intersects(m, q)
    @test intersects(q, m)
    @test intersects(m, m)

    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    b = Ball(P2(0, 0), T(1))
    @test intersects(t, t)
    @test intersects(b, b)
    @test intersects(t, b)
    @test intersects(b, t)

    t = Triangle(P2(1, 0), P2(2, 0), P2(1, 1))
    b = Ball(P2(0, 0), T(1))
    @test intersects(t, t)
    @test intersects(b, b)
    @test intersects(t, b)
    @test intersects(b, t)

    t = Triangle(P2(1, 0), P2(2, 0), P2(1, 1))
    b = Ball(P2(-0.01, 0), T(1))
    @test intersects(t, t)
    @test intersects(b, b)
    @test !intersects(t, b)
    @test !intersects(b, t)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/250
    t1 = Triangle(P3(0, 0, 0), P3(2, 0, 0), P3(1, 2, 0))
    t2 = Triangle(P3(1, 0, 0), P3(3, 0, 0), P3(2, 2, 0))
    t3 = Triangle(P3(3, 0, 0), P3(5, 0, 0), P3(4, 2, 0))
    @test intersects(t1, t2)
    @test intersects(t2, t3)
    @test !intersects(t1, t3)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/639
    r = Ray(P2(0.41169768366272996, 0.8990554132423699), V2(0.47249211625247445, 0.2523149692768657))
    b = Box(P2(1.0, 1.0), P2(5.0, 2.0))
    @test intersects(r, b)
    @test intersects(b, r)

    t = Triangle(P3(0, 0, 0), P3(2, 0, 0), P3(1, 2, 0))
    r1 = Ray(P3(1, 1, 1), V3(0, 0, -1))
    r2 = Ray(P3(1, 1, 1), V3(0, 0, 1))
    @test intersects(r1, t)
    @test intersects(t, r1)
    @test !intersects(r2, t)
    @test !intersects(t, r2)

    r = Ray(P2(0, 0), V2(1, 0))
    s1 = Sphere(P2(3, 0), T(1))
    s2 = Sphere(P2(0, 3), T(1))
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

    r = Ray(P2(0, 0), V2(1, 0))
    s = Sphere(P2(floatmax(Float32) / 2, 0), 1)
    @test intersects(r, s)

    r = Ray(P3(0, 0, 0), V3(1, 0, 0))
    s1 = Sphere(P3(5, 0, 1 - eps(T(1))), T(1))
    s2 = Sphere(P3(5, 0, 1 + eps(T(1))), T(1))
    @test intersects(r, s1)
    @test !intersects(r, s2)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/635
    q1 = Quadrangle(P3(4.0, 4.0, 0.0), P3(3.0, 3.0, 2.0), P3(3.0, 1.0, 2.0), P3(4.0, 0.0, 0.0))
    q2 = Quadrangle(P3(3.6, 3.0, 1.0), P3(5.6, 3.0, 1.0), P3(5.6, 1.0, 1.0), P3(3.6, 1.0, 1.0))
    q3 = Quadrangle(P3(3.6, 1.0, 1.0), P3(5.6, 1.0, 1.0), P3(5.6, -1.0, 1.0), P3(3.6, -1.0, 1.0))
    q4 = Quadrangle(P3(2.1, 1.0, 1.0), P3(4.1, 1.0, 1.0), P3(4.1, -1.0, 1.0), P3(2.1, -1.0, 1.0))
    @test !intersects(q1, q2)
    @test !intersects(q1, q3)
    @test intersects(q1, q1)
    @test intersects(q1, q4)

    h1 = Tetrahedron(P3(1, 1, 0), P3(4, 4, 0), P3(2.5, 2.5, 1.5), P3(1, 3, 2))
    h2 = Tetrahedron(P3(-1.0, 2.0, 1.0), P3(2.0, 1.0, 1.0), P3(-1.0, 4.0, 0.0), P3(0.5, 2.5, 1.5))
    h3 = Tetrahedron(P3(-1.3, 2.0, 1.0), P3(1.7, 1.0, 1.0), P3(-1.3, 4.0, 0.0), P3(0.2, 2.5, 1.5))
    @test intersects(h1, h2)
    @test !intersects(h1, h3)

    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly1 = PolyArea(outer)
    poly2 = PolyArea([outer, hole1, hole2])
    ball1 = Ball(P2(0.5, 0.5), T(0.05))
    ball2 = Ball(P2(0.3, 0.3), T(0.05))
    ball3 = Ball(P2(0.7, 0.3), T(0.05))
    ball4 = Ball(P2(0.3, 0.3), T(0.15))
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
    mesh1 = discretize(poly1, Dehn1899())
    mesh2 = discretize(poly2, Dehn1899())
    @test intersects(mesh1, mesh1)
    @test intersects(mesh2, mesh2)
    @test intersects(mesh1, mesh2)
    @test intersects(mesh2, mesh1)

    point = P2(0.5, 0.5)
    ball = Ball(P2(0, 0), T(1))
    @test intersects(point, ball)
    @test intersects(ball, point)
    @test intersects(point, point)
    @test !intersects(point, point + V2(1, 1))

    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    box = Box(P2(0, 0), P2(2, 2))
    @test intersects(poly, box)

    b1 = Box(P2(0, 0), P2(2, 2))
    b2 = Box(P2(2, 0), P2(4, 2))
    p1 = P2(1, 1)
    p2 = P2(3, 1)
    m = Multi([b1, b2])
    @test intersects(p1, b1)
    @test !intersects(p2, b1)
    @test intersects(p2, b2)
    @test !intersects(p1, b2)
    @test intersects(m, p1)
    @test intersects(p1, m)
    @test intersects(m, p2)
    @test intersects(p2, m)

    s1 = Segment(P2(0, 0), P2(4, 4))
    s2 = Segment(P2(4, 0), P2(0, 4))
    s3 = Segment(P2(2, 0), P2(4, 2))
    @test intersects(s1, s2)
    @test intersects(s2, s3)
    @test !intersects(s1, s3)

    s1 = Segment(P2(4, 0), P2(0, 4))
    s2 = Segment(P2(4, 0), P2(8, 4))
    s3 = Segment(P2(0, 8), P2(8, 8))
    r1 = Rope(P2[(0, 0), (4, 4), (8, 0)])
    r2 = Ring(P2[(0, 2), (4, 6), (8, 2)])
    @test intersects(s1, r1)
    @test intersects(s2, r1)
    @test !intersects(s3, r1)
    @test intersects(s1, r2)
    @test intersects(s2, r2)
    @test !intersects(s3, r2)
    @test intersects(r1, r2)

    r1 = Rope(P2[(0, 0), (2, 2), (4, 0)])
    r2 = Rope(P2[(3, 0), (5, 2), (7, 0)])
    r3 = Rope(P2[(6, 0), (8, 2), (10, 0)])
    @test intersects(r1, r2)
    @test intersects(r2, r3)
    @test !intersects(r1, r3)

    r1 = Ring(P2[(0, 0), (2, 2), (4, 0)])
    r2 = Ring(P2[(3, 0), (5, 2), (7, 0)])
    r3 = Ring(P2[(6, 0), (8, 2), (10, 0)])
    @test intersects(r1, r2)
    @test intersects(r2, r3)
    @test !intersects(r1, r3)

    t = Triangle(P2(3, 1), P2(7, 5), P2(11, 1))
    q = Quadrangle(P2(2, 0), P2(2, 7), P2(12, 7), P2(12, 0))
    b = Box(P2(2, 0), P2(12, 7))
    s1 = Segment(P2(5, 2), P2(9, 2))
    s2 = Segment(P2(0, 3), P2(5, 3))
    s3 = Segment(P2(4, 4), P2(10, 4))
    s4 = Segment(P2(1, 6), P2(13, 6))
    s5 = Segment(P2(0, 9), P2(14, 9))
    r1 = Ring(P2[(1, 2), (7, 8), (13, 2)])
    r2 = Rope(P2[(1, 2), (7, 8), (13, 2)])
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
    b1 = Box(P2(0, 0), P2(3, 3))
    b2 = Box(P2(2, 2), P2(5, 5))
    @test intersects(b1, b2)
    @test intersects(b2, b1)
    @test @elapsed(intersects(b1, b2)) < 5e-5
    @test @allocated(intersects(b1, b2)) < 100

    # partial application
    points = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    poly = PolyArea(points)
    box = Box(P2(0, 0), P2(2, 2))
    @test intersects(box)(poly)
    @test all(intersects(box), points)

    # method ambiguities
    point = P2(3, 1)
    ring = Ring(P2[(0, 0), (2, 2), (4, 0)])
    rope = Rope(P2[(2, 0), (4, 2), (6, 0)])
    seg = Segment(P2(0, 1), P2(6, 1))
    multi = Multi([ring])
    @test intersects(point, ring)
    @test intersects(point, rope)
    @test intersects(point, seg)
    @test intersects(point, multi)
    @test intersects(ring, multi)
    @test intersects(rope, multi)
    @test intersects(seg, multi)
  end

  @testset "iscollinear" begin
    @test iscollinear(P2(0, 0), P2(1, 1), P2(2, 2))
  end

  @testset "iscoplanar" begin
    @test iscoplanar(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0))
  end
end
