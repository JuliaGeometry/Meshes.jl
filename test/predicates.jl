@testset "Predicates" begin
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

    s1 = Segment(P2(0, 0), P2(1, 1))
    s2 = Segment(P2(0.5, 0.5), P2(1, 1))
    s3 = Segment(P2(0, 0), P2(0.5, 0.5))
    @test s2 ⊆ s1
    @test s3 ⊆ s1

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

    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly1 = PolyArea(outer)
    poly2 = PolyArea(outer, [hole1, hole2])
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
end
