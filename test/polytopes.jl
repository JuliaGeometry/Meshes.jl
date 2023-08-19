@testset "Polytopes" begin
  @testset "Segments" begin
    @test paramdim(Segment) == 1
    @test nvertices(Segment) == 2
    @test isperiodic(Segment) == (false,)
    @test isparametrized(Segment)

    s = Segment(P1(1.0), P1(2.0))
    @test vertex(s, 1) == P1(1.0)
    @test vertex(s, 2) == P1(2.0)
    @test all(P1(x) ∈ s for x in 1:0.01:2)
    @test all(P1(x) ∉ s for x in [-1.0, 0.0, 0.99, 2.1, 5.0, 10.0])
    @test s ≈ s
    @test !(s ≈ Segment(P1(2.0), P1(1.0)))
    @test !(s ≈ Segment(P1(-1.0), P1(2.0)))

    s = Segment(P2(0, 0), P2(1, 1))
    @test isconvex(s)
    @test isperiodic(s) == (false,)
    @test isparametrized(s)
    @test minimum(s) == P2(0, 0)
    @test maximum(s) == P2(1, 1)
    @test extrema(s) == (P2(0, 0), P2(1, 1))
    @test isapprox(length(s), sqrt(T(2)))
    @test s(T(0)) == P2(0, 0)
    @test s(T(1)) == P2(1, 1)
    @test all(P2(x, x) ∈ s for x in 0:0.01:1)
    @test all(p ∉ s for p in [P2(-0.1, -0.1), P2(1.1, 1.1), P2(0.5, 0.49), P2(1, 2)])
    @test_throws DomainError(T(1.2), "s(t) is not defined for t outside [0, 1].") s(T(1.2))
    @test_throws DomainError(T(-0.5), "s(t) is not defined for t outside [0, 1].") s(T(-0.5))
    @test s ≈ s
    @test !(s ≈ Segment(P2(1, 1), P2(0, 0)))
    @test !(s ≈ Segment(P2(1, 2), P2(0, 0)))

    s = Segment(P3(0, 0, 0), P3(1, 1, 1))
    @test all(P3(x, x, x) ∈ s for x in 0:0.01:1)
    @test all(p ∉ s for p in [P3(-0.1, -0.1, -0.1), P3(1.1, 1.1, 1.1)])
    @test all(p ∉ s for p in [P3(0.5, 0.5, 0.49), P3(1, 1, 2)])
    @test s ≈ s
    @test !(s ≈ Segment(P3(1, 1, 1), P3(0, 0, 0)))
    @test !(s ≈ Segment(P3(1, 1, 1), P3(0, 1, 0)))

    s = Segment(Point(1.0, 1.0, 1.0, 1.0), Point(2.0, 2.0, 2.0, 2.0))
    @test all(Point(x, x, x, x) ∈ s for x in 1:0.01:2)
    @test all(p ∉ s for p in [Point(0.99, 0.99, 0.99, 0.99), Point(2.1, 2.1, 2.1, 2.1)])
    @test all(p ∉ s for p in [Point(1.5, 1.5, 1.5, 1.49), Point(1, 1, 2, 1.0)])
    @test s ≈ s
    @test !(s ≈ Segment(Point(2, 2, 2, 2), Point(1, 1, 1, 1)))
    @test !(s ≈ Segment(Point(1, 1, 2, 1), Point(0, 0, 0, 0)))

    s = Segment(P3(0, 0, 0), P3(1, 1, 1))
    @test boundary(s) == PointSet([P3(0, 0, 0), P3(1, 1, 1)])
    @test perimeter(s) == zero(T)
    @test center(s) == P3(0.5, 0.5, 0.5)
    @test coordtype(center(s)) == T

    s = rand(Segment{2,T})
    @test s isa Segment
    @test embeddim(s) == 2
    @test coordtype(s) === T
    s = rand(Segment{3,T})
    @test s isa Segment
    @test embeddim(s) == 3
    @test coordtype(s) === T
  end

  @testset "Ropes/Rings" begin
    c1 = Rope(P2[(1, 1), (2, 2)])
    c2 = Rope(P2(1, 1), P2(2, 2))
    c3 = Rope(T.((1, 1.0)), T.((2.0, 2.0)))
    @test c1 == c2 == c3
    c1 = Ring(P2[(1, 1), (2, 2)])
    c2 = Ring(P2(1, 1), P2(2, 2))
    c3 = Ring(T.((1, 1.0)), T.((2.0, 2.0)))
    @test c1 == c2 == c3

    c = Rope(P2[(1, 1), (2, 2)])
    @test vertex(c, 1) == P2(1, 1)
    @test vertex(c, 2) == P2(2, 2)
    c = Ring(P2[(1, 1), (2, 2)])
    @test vertex(c, 0) == P2(2, 2)
    @test vertex(c, 1) == P2(1, 1)
    @test vertex(c, 2) == P2(2, 2)
    @test vertex(c, 3) == P2(1, 1)
    @test vertex(c, 4) == P2(2, 2)

    c = Rope(P2[(1, 1), (2, 2), (3, 3)])
    @test collect(segments(c)) == [Segment(P2(1, 1), P2(2, 2)), Segment(P2(2, 2), P2(3, 3))]
    c = Ring(P2[(1, 1), (2, 2), (3, 3)])
    @test collect(segments(c)) ==
          [Segment(P2(1, 1), P2(2, 2)), Segment(P2(2, 2), P2(3, 3)), Segment(P2(3, 3), P2(1, 1))]

    c = Rope(P2[(1, 1), (2, 2), (2, 2), (3, 3)])
    @test unique(c) == Rope(P2[(1, 1), (2, 2), (3, 3)])
    @test c == Rope(P2[(1, 1), (2, 2), (2, 2), (3, 3)])
    unique!(c)
    @test c == Rope(P2[(1, 1), (2, 2), (3, 3)])

    c = Rope(P2[(1, 1), (2, 2), (3, 3)])
    @test close(c) == Ring(P2[(1, 1), (2, 2), (3, 3)])
    c = Ring(P2[(1, 1), (2, 2), (3, 3)])
    @test open(c) == Rope(P2[(1, 1), (2, 2), (3, 3)])

    c = Rope(P2[(1, 1), (2, 2), (3, 3)])
    reverse!(c)
    @test c == Rope(P2[(3, 3), (2, 2), (1, 1)])
    c = Rope(P2[(1, 1), (2, 2), (3, 3)])
    @test reverse(c) == Rope(P2[(3, 3), (2, 2), (1, 1)])

    c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test angles(c) ≈ [-π / 2, -π / 2, -π / 2, -π / 2]
    c = Rope(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test angles(c) ≈ [-π / 2, -π / 2]
    c = Ring(P2[(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)])
    @test angles(c) ≈ [-atan(2), -π / 2, +π / 2, -π / 2, -π / 2, -(π - atan(2))]
    @test innerangles(c) ≈ [atan(2), π / 2, 3π / 2, π / 2, π / 2, π - atan(2)]

    c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test windingnumber(P2(0.5, 0.5), c) ≈ 1
    @test windingnumber(P2(0.5, 0.5), reverse(c)) ≈ -1
    c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0), (1, 0), (1, 1), (0, 1)])
    @test windingnumber(P2(0.5, 0.5), c) ≈ 2
    @test windingnumber(P2(0.5, 0.5), reverse(c)) ≈ -2

    c1 = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    c2 = Ring(vertices(c1))
    @test c1 == c2

    c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test centroid(c) == P2(0.5, 0.5)

    c = Rope(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test boundary(c) == PointSet(P2[(0, 0), (0, 1)])
    c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test isnothing(boundary(c))

    # should not repeat the first vertex manually
    @test_throws ArgumentError Ring(P2[(0, 0), (0, 0)])
    @test_throws ArgumentError Ring(P2[(0, 0), (1, 0), (1, 1), (0, 0)])

    # degenerate rings with 1 or 2 vertices are allowed
    r = Ring(P2[(0, 0)])
    @test isclosed(r)
    @test nvertices(r) == 1
    @test collect(segments(r)) == [Segment(P2(0, 0), P2(0, 0))]
    r = Ring(P2[(0, 0), (1, 1)])
    @test isclosed(r)
    @test nvertices(r) == 2
    @test collect(segments(r)) == [Segment(P2(0, 0), P2(1, 1)), Segment(P2(1, 1), P2(0, 0))]

    # orientation of 3D rings in X-Y plane
    r1 = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    r2 = Ring(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
    @test orientation(r1) == orientation(r2)

    p1 = P2(1, 1)
    p2 = P2(3, 1)
    p3 = P2(1, 0)
    p4 = P2(3, 0)
    pts = P2[(0, 0), (2, 2), (4, 0)]
    r = Ring(pts)
    @test p1 ∈ r
    @test p2 ∈ r
    @test p3 ∈ r
    @test p4 ∈ r
    r = Rope(pts)
    @test p1 ∈ r
    @test p2 ∈ r
    @test p3 ∉ r
    @test p4 ∉ r

    # approximately equal vertices
    pts = P2[
      (-48.04448403189499, -18.326530800015174)
      (-48.044478457836675, -18.326503670869467)
      (-48.04447845783733, -18.326503670869915)
      (-48.04447835073269, -18.326503149587666)
      (-48.044468448930644, -18.326490894176693)
      (-48.04447208741723, -18.326486301018672)
      (-48.044459173572015, -18.32646700775326)
      (-48.04445616736389, -18.326461847186216)
      (-48.044459897846174, -18.326466190774774)
      (-48.044462696066695, -18.32646303439271)
      (-48.044473299571635, -18.326478565399572)
      (-48.044473299571635, -18.326478565399565)
      (-48.044484052460334, -18.326494315209573)
      (-48.04449288424675, -18.326504598503668)
      (-48.044492356262886, -18.32650647783081)
      (-48.0444943180541, -18.326509351276243)
      (-48.044492458690776, -18.32651322842786)
      (-48.04450917793127, -18.326524641668517)
      (-48.044501408820125, -18.326551273900744)
    ]
    r1 = Rope(pts)
    r2 = Ring(pts)
    ur1 = unique(r1)
    ur2 = unique(r2)
    @test nvertices(ur1) < nvertices(r1)
    @test nvertices(ur2) < nvertices(r2)
    if T === Float32
      @test nvertices(ur1) == 10
      @test nvertices(ur2) == 10
    else
      @test nvertices(ur1) == 17
      @test nvertices(ur2) == 17
    end

    r = rand(Rope{2,T})
    @test r isa Rope
    @test embeddim(r) == 2
    @test coordtype(r) === T
    r = rand(Rope{3,T})
    @test r isa Rope
    @test embeddim(r) == 3
    @test coordtype(r) === T

    r = rand(Ring{2,T})
    @test r isa Ring
    @test embeddim(r) == 2
    @test coordtype(r) === T
    r = rand(Ring{3,T})
    @test r isa Ring
    @test embeddim(r) == 3
    @test coordtype(r) === T

    # issimple benchmark
    r = Sphere(P2(0, 0), T(1)) |> pointify |> Ring
    @test issimple(r)
    @test @elapsed(issimple(r)) < 0.02
    @test @allocated(issimple(r)) < 950000

    # innerangles in 3D is obtained via projection
    r1 = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    r2 = Ring(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
    @test innerangles(r1) ≈ innerangles(r2)
  end

  @testset "Ngons" begin
    @test paramdim(Ngon) == 2
    NGONS = [Triangle, Quadrangle, Pentagon, Hexagon, Heptagon, Octagon, Nonagon, Decagon]
    NVERT = 3:10
    for (i, NGON) in enumerate(NGONS)
      @test paramdim(NGON) == 2
      @test nvertices(NGON) == NVERT[i]

      n = rand(NGON{2,T})
      @test n isa NGON
      @test embeddim(n) == 2
      @test coordtype(n) === T
      n = rand(NGON{3,T})
      @test n isa NGON
      @test embeddim(n) == 3
      @test coordtype(n) === T
    end

    # ---------
    # TRIANGLE
    # ---------

    @test issimplex(Triangle)
    @test isparametrized(Triangle)

    # Triangle in 2D space
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test issimplex(t)
    @test isconvex(t)
    @test isparametrized(t)
    @test vertex(t, 1) == P2(0, 0)
    @test vertex(t, 2) == P2(1, 0)
    @test vertex(t, 3) == P2(0, 1)
    @test signarea(t) == T(0.5)
    @test area(t) == T(0.5)
    t = Triangle(P2(0, 0), P2(0, 1), P2(1, 0))
    @test signarea(t) == T(-0.5)
    @test area(t) == T(0.5)
    t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
    for p in P2[(0, 0), (1, 0), (1, 1), (0.5, 0.0), (1.0, 0.5), (0.5, 0.5)]
      @test p ∈ t
    end
    for p in P2[(-1, 0), (0, -1), (0.5, 1.0)]
      @test p ∉ t
    end
    t = Triangle(P2(0.4, 0.4), P2(0.6, 0.4), P2(0.8, 0.4))
    @test P2(0.2, 0.4) ∉ t
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test t(T(0.0), T(0.0)) == P2(0, 0)
    @test t(T(1.0), T(0.0)) == P2(1, 0)
    @test t(T(0.0), T(1.0)) == P2(0, 1)
    @test t(T(0.5), T(0.5)) == P2(0.5, 0.5)
    @test_throws DomainError((T(-0.5), T(0.0)), "invalid barycentric coordinates for triangle.") t(T(-0.5), T(0.0))
    @test_throws DomainError((T(1), T(1)), "invalid barycentric coordinates for triangle.") t(T(1), T(1))
    @test !hasholes(t)
    @test unique(t) == t
    @test boundary(t) == first(rings(t))
    @test rings(t) == [Ring(P2(0, 0), P2(1, 0), P2(0, 1))]
    @test (t |> Bridge() |> boundary) == boundary(t)

    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test perimeter(t) ≈ T(1 + 1 + √2)

    # https://github.com/JuliaGeometry/Meshes.jl/issues/333
    t = Triangle((0.0f0, 0.0f0), (1.0f0, 0.0f0), (0.5f0, 1.0f0))
    @test Point(0.5f0, 0.5f0) ∈ t
    @test Point(0.5e0, 0.5e0) ∈ t

    # test orientation
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test orientation(t) == :CCW
    t = Triangle(P2(0, 0), P2(0, 1), P2(1, 0))
    @test orientation(t) == :CW

    # test angles
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test all(isapprox.(rad2deg.(angles(t)), T[-90, -45, -45], atol=8 * eps(T)))
    @test all(isapprox.(rad2deg.(innerangles(t)), T[90, 45, 45], atol=8 * eps(T)))

    # Triangle in 3D space
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    @test area(t) == T(0.5)
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 1))
    @test area(t) > T(0.7)
    for p in P3[(0, 0, 0), (1, 0, 0), (0, 1, 1), (0, 0.2, 0.2)]
      @test p ∈ t
    end
    for p in P3[(-1, 0, 0), (1, 2, 0), (0, 1, 2)]
      @test p ∉ t
    end
    t = Triangle(P3(0, 0, 0), P3(0, 1, 0), P3(0, 0, 1))
    @test t(T(0.0), T(0.0)) == P3(0, 0, 0)
    @test t(T(1.0), T(0.0)) == P3(0, 1, 0)
    @test t(T(0.0), T(1.0)) == P3(0, 0, 1)
    @test t(T(0.5), T(0.5)) == P3(0, 0.5, 0.5)
    @test_throws DomainError((T(-0.5), T(0.0)), "invalid barycentric coordinates for triangle.") t(T(-0.5), T(0.0))
    @test_throws DomainError((T(1), T(1)), "invalid barycentric coordinates for triangle.") t(T(1), T(1))
    @test isapprox(normal(t), Vec(1, 0, 0))
    t = Triangle(P3(0, 0, 0), P3(2, 0, 0), P3(0, 2, 2))
    @test isapprox(normal(t), Vec(0, -1 / sqrt(2), 1 / sqrt(2)))

    # test convexity of Triangle
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    @test isconvex(t)
    t = Triangle(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0))
    @test isconvex(t)

    # -----------
    # QUADRANGLE
    # -----------

    @test isperiodic(Quadrangle) == (false, false)
    @test isparametrized(Quadrangle)

    # test periodicity of Quadrangle
    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test isperiodic(q) == (false, false)
    @test isparametrized(q)

    # Quadrangle in 2D space
    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test vertex(q, 1) == P2(0, 0)
    @test vertex(q, 2) == P2(1, 0)
    @test vertex(q, 3) == P2(1, 1)
    @test vertex(q, 4) == P2(0, 1)
    @test area(q) == T(1)
    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1.5, 1.0), P2(0.5, 1.0))
    @test area(q) == T(1)
    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1.5, 1.0), P2(0.5, 1.0))
    for p in P2[(0, 0), (1, 0), (1.5, 1.0), (0.5, 1.0), (0.5, 0.5)]
      @test p ∈ q
    end
    for p in P2[(0, 1), (1.5, 0.0)]
      @test p ∉ q
    end
    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test !hasholes(q)
    @test unique(q) == q
    @test boundary(q) == first(rings(q))
    @test rings(q) == [Ring(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))]
    @test (q |> Bridge() |> boundary) == boundary(q)
    @test q(T(0), T(0)) == P2(0, 0)
    @test q(T(1), T(0)) == P2(1, 0)
    @test q(T(1), T(1)) == P2(1, 1)
    @test q(T(0), T(1)) == P2(0, 1)

    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test_throws DomainError((T(1.2), T(1.2)), "q(u, v) is not defined for u, v outside [0, 1]².") q(T(1.2), T(1.2))

    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test perimeter(q) ≈ T(4)

    # Quadrangle in 3D space
    q = Quadrangle(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0))
    @test area(q) == T(1)
    q = Quadrangle(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 1))
    @test area(q) > T(1)
    @test q(T(0), T(0)) == P3(0, 0, 0)
    @test q(T(1), T(0)) == P3(1, 0, 0)
    @test q(T(1), T(1)) == P3(1, 1, 0)
    @test q(T(0), T(1)) == P3(0, 1, 1)

    # isconvex in 2D
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

    # isconvex in 3D
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
  end

  @testset "PolyAreas" begin
    @test paramdim(PolyArea) == 2

    # equality and approximate equality
    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly = PolyArea(outer, [hole1, hole2])
    @test poly == poly
    @test poly ≈ poly

    # outer chain with 2 vertices is fixed by default
    poly = PolyArea(P2[(0, 0), (1, 0)])
    @test rings(poly) == [Ring(P2[(0, 0), (0.5, 0.0), (1, 0)])]

    # inner chain with 2 vertices is removed by default
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)], [P2[(1, 2), (2, 3)]])
    @test rings(poly) == [Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])]

    # orientation of chains is fixed by default
    poly = PolyArea(P2[(0, 0), (0, 1), (1, 1), (1, 0)])
    @test vertices(poly) == CircularVector(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    poly = PolyArea(P2[(0, 0), (0, 1), (1, 1), (1, 0)], fix=false)
    @test vertices(poly) == CircularVector(P2[(0, 0), (0, 1), (1, 1), (1, 0)])

    # test accessor methods
    poly = PolyArea(P2[(1, 2), (2, 3)], fix=false)
    @test vertices(poly) == CircularVector(P2[(1, 2), (2, 3)])
    poly = PolyArea(P2[(1, 2), (2, 3)], [P2[(1.1, 2.54), (1.4, 1.5)]], fix=false)
    @test vertices(poly) == CircularVector(P2[(1, 2), (2, 3), (1.1, 2.54), (1.4, 1.5)])

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output poly1
    fnames = ["poly$i.line" for i in 1:5]
    polys1 = [readpoly(T, joinpath(datadir, fname)) for fname in fnames]
    for poly in polys1
      @test !hasholes(poly)
      @test issimple(poly)
      @test boundary(poly) == first(rings(poly))
      @test nvertices(poly) == 30
      for algo in [WindingOrientation(), TriangleOrientation()]
        @test orientation(poly, algo) == :CCW
      end
      @test unique(poly) == poly
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output smooth1 --smooth 2
    fnames = ["smooth$i.line" for i in 1:5]
    polys2 = [readpoly(T, joinpath(datadir, fname)) for fname in fnames]
    for poly in polys2
      @test !hasholes(poly)
      @test issimple(poly)
      @test boundary(poly) == first(rings(poly))
      @test nvertices(poly) == 120
      for algo in [WindingOrientation(), TriangleOrientation()]
        @test orientation(poly, algo) == :CCW
      end
      @test unique(poly) == poly
    end

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output hole1 --holes 2
    fnames = ["hole$i.line" for i in 1:5]
    polys3 = [readpoly(T, joinpath(datadir, fname)) for fname in fnames]
    for poly in polys3
      rs = rings(poly)
      @test hasholes(poly)
      @test !issimple(poly)
      @test boundary(poly) == Multi(rs)
      @test nvertices(first(rs)) < 30
      @test all(nvertices.(rs[2:end]) .< 18)
      for algo in [WindingOrientation(), TriangleOrientation()]
        orients = orientation(poly, algo)
        @test orients[1] == :CCW
        @test all(orients[2:end] .== :CW)
      end
      @test unique(poly) == poly
    end

    # test bridges
    for poly in [polys1; polys2; polys3]
      b = poly |> Bridge()
      nb = nvertices(b)
      np = nvertices.(rings(poly))
      @test nb ≥ sum(np)
      # triangle orientation always works even
      # in the presence of self-intersections
      @test orientation(b, TriangleOrientation()) == :CCW
      # winding orientation is only suitable
      # for simple polygonal chains
      if issimple(b)
        @test orientation(b, WindingOrientation()) == :CCW
      end
    end

    # bridges between holes
    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly = PolyArea(outer, [hole1, hole2])
    @test vertices(poly) == P2[
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
    ]
    chain = poly |> Bridge()
    target = P2[
      (0.0, 0.0),
      (0.2, 0.2),
      (0.2, 0.4),
      (0.4, 0.4),
      (0.6, 0.4),
      (0.8, 0.4),
      (0.8, 0.2),
      (0.6, 0.2),
      (0.6, 0.4),
      (0.4, 0.4),
      (0.4, 0.2),
      (0.2, 0.2),
      (0.0, 0.0),
      (1.0, 0.0),
      (1.0, 1.0),
      (0.0, 1.0)
    ]
    @test vertices(chain) == target

    # test uniqueness
    points = P2[(1, 1), (2, 2), (2, 2), (3, 3)]
    poly = PolyArea(points)
    unique!(poly)
    @test first(rings(poly)) == Ring(P2[(1, 1), (2, 2), (3, 3)])

    # unique and bridges
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1), (0, 1)])
    cpoly = poly |> Repair{0}() |> Bridge()
    @test cpoly == PolyArea(P2[(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)])

    # approximately equal vertices
    poly = PolyArea(
      P2[
        (-48.04448403189499, -18.326530800015174)
        (-48.044478457836675, -18.326503670869467)
        (-48.04447845783733, -18.326503670869915)
        (-48.04447835073269, -18.326503149587666)
        (-48.044468448930644, -18.326490894176693)
        (-48.04447208741723, -18.326486301018672)
        (-48.044459173572015, -18.32646700775326)
        (-48.04445616736389, -18.326461847186216)
        (-48.044459897846174, -18.326466190774774)
        (-48.044462696066695, -18.32646303439271)
        (-48.044473299571635, -18.326478565399572)
        (-48.044473299571635, -18.326478565399565)
        (-48.044484052460334, -18.326494315209573)
        (-48.04449288424675, -18.326504598503668)
        (-48.044492356262886, -18.32650647783081)
        (-48.0444943180541, -18.326509351276243)
        (-48.044492458690776, -18.32651322842786)
        (-48.04450917793127, -18.326524641668517)
        (-48.044501408820125, -18.326551273900744)
      ]
    )
    upoly = unique(poly)
    @test nvertices(upoly) < nvertices(poly)
    if T === Float32
      @test nvertices(upoly) == 10
    else
      @test nvertices(upoly) == 17
    end

    # invalid inner
    outer = rand(P2, 10)
    v1, v2 = rand(V2, 2)
    inner = [Point(v1), Point(v1), Point(v2)]
    poly = PolyArea(outer, [inner])
    upoly = unique(poly)
    @test hasholes(poly)
    @test !hasholes(upoly)

    # centroid
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test centroid(poly) == P2(0.5, 0.5)

    # single vertex access
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    @test vertex(poly, 1) == P2(0, 0)
    @test vertex(poly, 4) == P2(0, 1)

    # point in polygonal area
    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly = PolyArea(outer, [hole1, hole2])
    @test all(p ∈ poly for p in outer)
    @test P2(0.5, 0.5) ∈ poly
    @test P2(0.2, 0.6) ∈ poly
    @test P2(1.5, 0.5) ∉ poly
    @test P2(-0.5, 0.5) ∉ poly
    @test P2(0.25, 0.25) ∉ poly
    @test P2(0.75, 0.25) ∉ poly
    @test P2(0.75, 0.75) ∈ poly

    # area
    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly = PolyArea(outer, [hole1, hole2])
    @test area(poly) ≈ T(0.92)

    # convexity
    outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
    hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
    hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
    poly1 = PolyArea(outer)
    poly2 = PolyArea(outer, [hole1, hole2])
    @test isconvex(poly1)
    @test !isconvex(poly2)
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0.5, 0.5), (0, 1)])
    @test !isconvex(poly)

    p = rand(PolyArea{2,T})
    @test p isa PolyArea
    @test embeddim(p) == 2
    @test coordtype(p) === T
    p = rand(PolyArea{3,T})
    @test p isa PolyArea
    @test embeddim(p) == 3
    @test coordtype(p) === T

    # should not repeat the first vertex manually
    @test_throws ArgumentError PolyArea(P2[(0, 0), (0, 0)])
    @test_throws ArgumentError PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 0)])
  end

  @testset "Polyhedra" begin
    @test paramdim(Tetrahedron) == 3
    @test nvertices(Tetrahedron) == 4

    t = Tetrahedron(P3(0, 0, 0), P3(1, 0, 0), P3(0, 1, 0), P3(0, 0, 1))
    @test vertex(t, 1) == P3(0, 0, 0)
    @test vertex(t, 2) == P3(1, 0, 0)
    @test vertex(t, 3) == P3(0, 1, 0)
    @test vertex(t, 4) == P3(0, 0, 1)
    @test issimplex(t)
    @test isconvex(t)
    @test measure(t) == T(1 / 6)
    m = boundary(t)
    n = normal.(m)
    @test m isa Mesh
    @test nvertices(m) == 4
    @test nelements(m) == 4
    @test n[1] == T[0, 0, -1]
    @test n[2] == T[0, -1, 0]
    @test n[3] == T[-1, 0, 0]
    @test all(>(0), n[4])
    @test t(T(0), T(0), T(0)) ≈ P3(0, 0, 0)
    @test t(T(1), T(0), T(0)) ≈ P3(1, 0, 0)
    @test t(T(0), T(1), T(0)) ≈ P3(0, 1, 0)
    @test t(T(0), T(0), T(1)) ≈ P3(0, 0, 1)
    @test_throws DomainError((T(1), T(1), T(1)), "invalid barycentric coordinates for tetrahedron.") t(T(1), T(1), T(1))

    t = rand(Tetrahedron{3,T})
    @test t isa Tetrahedron
    @test embeddim(t) == 3
    @test coordtype(t) === T

    @test paramdim(Hexahedron) == 3
    @test nvertices(Hexahedron) == 8
    @test isperiodic(Hexahedron) == (false, false, false)
    @test isparametrized(Hexahedron)

    h =
      Hexahedron(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0), P3(0, 0, 1), P3(1, 0, 1), P3(1, 1, 1), P3(0, 1, 1))
    @test vertex(h, 1) == P3(0, 0, 0)
    @test vertex(h, 8) == P3(0, 1, 1)
    @test isperiodic(h) == (false, false, false)
    @test isparametrized(h)
    @test h(T(0), T(0), T(0)) == P3(0, 0, 0)
    @test h(T(0), T(0), T(1)) == P3(0, 0, 1)
    @test h(T(0), T(1), T(0)) == P3(0, 1, 0)
    @test h(T(0), T(1), T(1)) == P3(0, 1, 1)
    @test h(T(1), T(0), T(0)) == P3(1, 0, 0)
    @test h(T(1), T(0), T(1)) == P3(1, 0, 1)
    @test h(T(1), T(1), T(0)) == P3(1, 1, 0)
    @test h(T(1), T(1), T(1)) == P3(1, 1, 1)

    h =
      Hexahedron(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0), P3(0, 0, 1), P3(1, 0, 1), P3(1, 1, 1), P3(0, 1, 1))
    @test volume(h) ≈ T(1 * 1 * 1)
    h =
      Hexahedron(P3(0, 0, 0), P3(2, 0, 0), P3(2, 2, 0), P3(0, 2, 0), P3(0, 0, 2), P3(2, 0, 2), P3(2, 2, 2), P3(0, 2, 2))
    @test volume(h) ≈ T(2 * 2 * 2)

    # volume formula of a frustum of a prism is V = 1/3*H*(S₁+S₂+sqrt(S₁*S₂))
    # here we build a hexahedron which is a frustum of a prism with
    # bottom area S₁= 4, top area S₂= 1, height H = 2
    h =
      Hexahedron(P3(0, 0, 0), P3(2, 0, 0), P3(2, 2, 0), P3(0, 2, 0), P3(0, 0, 2), P3(1, 0, 2), P3(1, 1, 2), P3(0, 1, 2))
    @test volume(h) ≈ T(1 / 3 * 2 * (1 + 4 + sqrt(1 * 4)))

    h =
      Hexahedron(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0), P3(0, 0, 1), P3(1, 0, 1), P3(1, 1, 1), P3(0, 1, 1))
    m = boundary(h)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    h = rand(Hexahedron{3,T})
    @test h isa Hexahedron
    @test embeddim(h) == 3
    @test coordtype(h) === T

    @test paramdim(Pyramid) == 3
    @test nvertices(Pyramid) == 5

    p = Pyramid(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0), P3(0, 0, 1))
    @test volume(p) ≈ T(1 / 3)
    m = boundary(p)
    @test m isa Mesh
    @test nelements(m) == 5
    @test m[1] isa Quadrangle
    @test m[2] isa Triangle
    @test m[3] isa Triangle
    @test m[4] isa Triangle
    @test m[5] isa Triangle

    p = rand(Pyramid{3,T})
    @test p isa Pyramid
    @test embeddim(p) == 3
    @test coordtype(p) === T
  end
end
