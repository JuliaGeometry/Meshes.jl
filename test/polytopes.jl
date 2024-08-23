@testitem "Polytopes" begin
  @testitem "Segments" begin
    @test paramdim(Segment) == 1
    @test nvertices(Segment) == 2

    s = Segment(cart(1.0), cart(2.0))
    @test crs(s) <: Cartesian{NoDatum}
    @test Meshes.lentype(s) == ℳ
    @test vertex(s, 1) == cart(1.0)
    @test vertex(s, 2) == cart(2.0)
    @test all(cart(x) ∈ s for x in 1:0.01:2)
    @test all(cart(x) ∉ s for x in [-1.0, 0.0, 0.99, 2.1, 5.0, 10.0])
    @test s ≈ s
    @test !(s ≈ Segment(cart(2.0), cart(1.0)))
    @test !(s ≈ Segment(cart(-1.0), cart(2.0)))
    @test reverse(s) == Segment(cart(2.0), cart(1.0))

    s = Segment(cart(0, 0), cart(1, 1))
    @test minimum(s) == cart(0, 0)
    @test maximum(s) == cart(1, 1)
    @test extrema(s) == (cart(0, 0), cart(1, 1))
    @test isapprox(length(s), sqrt(T(2)) * u"m")
    @test s(T(0)) == cart(0, 0)
    @test s(T(1)) == cart(1, 1)
    @test all(cart(x, x) ∈ s for x in 0:0.01:1)
    @test all(p ∉ s for p in [cart(-0.1, -0.1), cart(1.1, 1.1), cart(0.5, 0.49), cart(1, 2)])
    @test_throws DomainError(T(1.2), "s(t) is not defined for t outside [0, 1].") s(T(1.2))
    @test_throws DomainError(T(-0.5), "s(t) is not defined for t outside [0, 1].") s(T(-0.5))
    @test s ≈ s
    @test !(s ≈ Segment(cart(1, 1), cart(0, 0)))
    @test !(s ≈ Segment(cart(1, 2), cart(0, 0)))
    @test reverse(s) == Segment(cart(1, 1), cart(0, 0))

    s = Segment(cart(0, 0, 0), cart(1, 1, 1))
    @test all(cart(x, x, x) ∈ s for x in 0:0.01:1)
    @test all(p ∉ s for p in [cart(-0.1, -0.1, -0.1), cart(1.1, 1.1, 1.1)])
    @test all(p ∉ s for p in [cart(0.5, 0.5, 0.49), cart(1, 1, 2)])
    @test s ≈ s
    @test !(s ≈ Segment(cart(1, 1, 1), cart(0, 0, 0)))
    @test !(s ≈ Segment(cart(1, 1, 1), cart(0, 1, 0)))
    @test reverse(s) == Segment(cart(1, 1, 1), cart(0, 0, 0))

    s = Segment(cart(0, 0), cart(1, 1))
    equaltest(s)
    isapproxtest(s)

    s = Segment(Point(1.0, 1.0, 1.0, 1.0), Point(2.0, 2.0, 2.0, 2.0))
    @test all(Point(x, x, x, x) ∈ s for x in 1:0.01:2)
    @test all(p ∉ s for p in [Point(0.99, 0.99, 0.99, 0.99), Point(2.1, 2.1, 2.1, 2.1)])
    @test all(p ∉ s for p in [Point(1.5, 1.5, 1.5, 1.49), Point(1, 1, 2, 1.0)])
    @test s ≈ s
    @test !(s ≈ Segment(Point(2, 2, 2, 2), Point(1, 1, 1, 1)))
    @test !(s ≈ Segment(Point(1, 1, 2, 1), Point(0, 0, 0, 0)))

    s = Segment(cart(0, 0, 0), cart(1, 1, 1))
    @test boundary(s) == Multi([cart(0, 0, 0), cart(1, 1, 1)])
    @test perimeter(s) == zero(T) * u"m"
    @test centroid(s) == cart(0.5, 0.5, 0.5)
    @test Meshes.lentype(centroid(s)) == ℳ

    # unitful coordinates
    x1 = T(0)u"m"
    x2 = T(1)u"m"
    s = Segment(Point(x1, x1, x1), Point(x2, x2, x2))
    @test boundary(s) == Multi([Point(x1, x1, x1), Point(x2, x2, x2)])
    @test perimeter(s) == 0u"m"
    xm = T(0.5)u"m"
    @test centroid(s) == Point(xm, xm, xm)
    @test Meshes.lentype(centroid(s)) == typeof(xm)

    # CRS propagation
    s = Segment(merc(0, 0), merc(1, 1))
    @test crs(s(T(0))) === crs(s)

    s = Segment(cart(0, 0), cart(1, 1))
    @test sprint(show, s) == "Segment((x: 0.0 m, y: 0.0 m), (x: 1.0 m, y: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), s) == """
      Segment
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
      └─ Point(x: 1.0f0 m, y: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), s) == """
      Segment
      ├─ Point(x: 0.0 m, y: 0.0 m)
      └─ Point(x: 1.0 m, y: 1.0 m)"""
    end
  end

  @testitem "Ropes/Rings" begin
    c1 = Rope(cart.([(1, 1), (2, 2)]))
    c2 = Rope(cart(1, 1), cart(2, 2))
    c3 = Rope(T.((1, 1.0)), T.((2.0, 2.0)))
    @test c1 == c2 == c3
    c1 = Ring(cart.([(1, 1), (2, 2)]))
    c2 = Ring(cart(1, 1), cart(2, 2))
    c3 = Ring(T.((1, 1.0)), T.((2.0, 2.0)))
    @test c1 == c2 == c3

    c = Rope(cart(0, 0), cart(1, 0), cart(0, 1))
    equaltest(c)
    isapproxtest(c)

    c = Ring(cart(0, 0), cart(1, 0), cart(0, 1))
    equaltest(c)
    isapproxtest(c)

    # circular equality
    c1 = Ring(cart.([(1, 1), (2, 2), (3, 3)]))
    c2 = Ring(cart.([(2, 2), (3, 3), (1, 1)]))
    c3 = Ring(cart.([(3, 3), (1, 1), (2, 2)]))
    @test c1 ≗ c2 ≗ c3

    c = Rope(cart.([(1, 1), (2, 2)]))
    @test crs(c) <: Cartesian{NoDatum}
    @test Meshes.lentype(c) == ℳ
    @test vertex(c, 1) == cart(1, 1)
    @test vertex(c, 2) == cart(2, 2)
    c = Ring(cart.([(1, 1), (2, 2)]))
    @test crs(c) <: Cartesian{NoDatum}
    @test Meshes.lentype(c) == ℳ
    @test vertex(c, 0) == cart(2, 2)
    @test vertex(c, 1) == cart(1, 1)
    @test vertex(c, 2) == cart(2, 2)
    @test vertex(c, 3) == cart(1, 1)
    @test vertex(c, 4) == cart(2, 2)

    c = Rope(cart.([(1, 1), (2, 2), (3, 3)]))
    @test collect(segments(c)) == [Segment(cart(1, 1), cart(2, 2)), Segment(cart(2, 2), cart(3, 3))]
    c = Ring(cart.([(1, 1), (2, 2), (3, 3)]))
    @test collect(segments(c)) ==
          [Segment(cart(1, 1), cart(2, 2)), Segment(cart(2, 2), cart(3, 3)), Segment(cart(3, 3), cart(1, 1))]

    c = Rope(cart.([(1, 1), (2, 2), (2, 2), (3, 3)]))
    @test unique(c) == Rope(cart.([(1, 1), (2, 2), (3, 3)]))
    @test c == Rope(cart.([(1, 1), (2, 2), (2, 2), (3, 3)]))
    unique!(c)
    @test c == Rope(cart.([(1, 1), (2, 2), (3, 3)]))

    c = Rope(cart.([(1, 1), (2, 2), (3, 3)]))
    @test close(c) == Ring(cart.([(1, 1), (2, 2), (3, 3)]))
    c = Ring(cart.([(1, 1), (2, 2), (3, 3)]))
    @test open(c) == Rope(cart.([(1, 1), (2, 2), (3, 3)]))

    c = Rope(cart.([(1, 1), (2, 2), (3, 3)]))
    reverse!(c)
    @test c == Rope(cart.([(3, 3), (2, 2), (1, 1)]))
    c = Rope(cart.([(1, 1), (2, 2), (3, 3)]))
    @test reverse(c) == Rope(cart.([(3, 3), (2, 2), (1, 1)]))

    c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test angles(c) ≈ [-π / 2, -π / 2, -π / 2, -π / 2]
    c = Rope(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test angles(c) ≈ [-π / 2, -π / 2]
    c = Ring(cart.([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]))
    @test angles(c) ≈ [-atan(2), -π / 2, +π / 2, -π / 2, -π / 2, -(π - atan(2))]
    @test innerangles(c) ≈ [atan(2), π / 2, 3π / 2, π / 2, π / 2, π - atan(2)]

    c1 = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    c2 = Ring(vertices(c1))
    @test c1 == c2

    c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test centroid(c) == cart(0.5, 0.5)

    c = Rope(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test boundary(c) == Multi(cart.([(0, 0), (0, 1)]))
    c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test isnothing(boundary(c))

    # should not repeat the first vertex manually
    @test_throws ArgumentError Ring(cart.([(0, 0), (0, 0)]))
    @test_throws ArgumentError Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 0)]))

    # degenerate rings with 1 or 2 vertices are allowed
    r = Ring(cart.([(0, 0)]))
    @test isclosed(r)
    @test nvertices(r) == 1
    @test collect(segments(r)) == [Segment(cart(0, 0), cart(0, 0))]
    r = Ring(cart.([(0, 0), (1, 1)]))
    @test isclosed(r)
    @test nvertices(r) == 2
    @test collect(segments(r)) == [Segment(cart(0, 0), cart(1, 1)), Segment(cart(1, 1), cart(0, 0))]

    p1 = cart(1, 1)
    p2 = cart(3, 1)
    p3 = cart(1, 0)
    p4 = cart(3, 0)
    pts = cart.([(0, 0), (2, 2), (4, 0)])
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
    pts =
      cart.(
        [
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

    # issimple benchmark
    r = Sphere(cart(0, 0), T(1)) |> pointify |> Ring
    @test issimple(r)
    @test @elapsed(issimple(r)) < 0.02
    @test @allocated(issimple(r)) < 950000

    # CRS propagation
    r = Ring(merc.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test crs(centroid(r)) === crs(r)

    ri = Ring(cart.([(1, 1), (2, 2), (3, 3)]))
    ro = Rope(cart.([(1, 1), (2, 2), (3, 3)]))
    @test sprint(show, ri) == "Ring((x: 1.0 m, y: 1.0 m), (x: 2.0 m, y: 2.0 m), (x: 3.0 m, y: 3.0 m))"
    @test sprint(show, ro) == "Rope((x: 1.0 m, y: 1.0 m), (x: 2.0 m, y: 2.0 m), (x: 3.0 m, y: 3.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), ri) == """
      Ring
      ├─ Point(x: 1.0f0 m, y: 1.0f0 m)
      ├─ Point(x: 2.0f0 m, y: 2.0f0 m)
      └─ Point(x: 3.0f0 m, y: 3.0f0 m)"""
      @test sprint(show, MIME("text/plain"), ro) == """
      Rope
      ├─ Point(x: 1.0f0 m, y: 1.0f0 m)
      ├─ Point(x: 2.0f0 m, y: 2.0f0 m)
      └─ Point(x: 3.0f0 m, y: 3.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), ri) == """
      Ring
      ├─ Point(x: 1.0 m, y: 1.0 m)
      ├─ Point(x: 2.0 m, y: 2.0 m)
      └─ Point(x: 3.0 m, y: 3.0 m)"""
      @test sprint(show, MIME("text/plain"), ro) == """
      Rope
      ├─ Point(x: 1.0 m, y: 1.0 m)
      ├─ Point(x: 2.0 m, y: 2.0 m)
      └─ Point(x: 3.0 m, y: 3.0 m)"""
    end
  end

  @testitem "Ngons" begin
    pts = (cart(0, 0), cart(1, 0), cart(0, 1))
    tups = (T.((0, 0)), T.((1, 0)), T.((0, 1)))
    @test paramdim(Ngon) == 2
    @test vertices(Ngon(pts)) == pts
    @test vertices(Ngon(pts...)) == pts
    @test vertices(Ngon(tups...)) == pts
    @test vertices(Ngon{3}(pts)) == pts
    @test vertices(Ngon{3}(pts...)) == pts
    @test vertices(Ngon{3}(tups...)) == pts

    NGONS = [Triangle, Quadrangle, Pentagon, Hexagon, Heptagon, Octagon, Nonagon, Decagon]
    NVERT = 3:10
    for (i, NGON) in enumerate(NGONS)
      @test paramdim(NGON) == 2
      @test nvertices(NGON) == NVERT[i]
    end

    # error: the number of vertices must be greater than or equal to 3
    @test_throws ArgumentError Ngon(cart(0, 0), cart(1, 1))
    @test_throws ArgumentError Ngon{2}(cart(0, 0), cart(1, 1))

    # ---------
    # TRIANGLE
    # ---------

    # Triangle in 2D space
    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test crs(t) <: Cartesian{NoDatum}
    @test Meshes.lentype(t) == ℳ
    @test vertex(t, 1) == cart(0, 0)
    @test vertex(t, 2) == cart(1, 0)
    @test vertex(t, 3) == cart(0, 1)
    @test area(t) == T(0.5) * u"m^2"
    t = Triangle(cart(0, 0), cart(0, 1), cart(1, 0))
    @test area(t) == T(0.5) * u"m^2"
    t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
    for p in cart.([(0, 0), (1, 0), (1, 1), (0.5, 0.0), (1.0, 0.5), (0.5, 0.5)])
      @test p ∈ t
    end
    for p in cart.([(-1, 0), (0, -1), (0.5, 1.0)])
      @test p ∉ t
    end
    t = Triangle(cart(0.4, 0.4), cart(0.6, 0.4), cart(0.8, 0.4))
    @test cart(0.2, 0.4) ∉ t
    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test t(T(0.0), T(0.0)) == cart(0, 0)
    @test t(T(1.0), T(0.0)) == cart(1, 0)
    @test t(T(0.0), T(1.0)) == cart(0, 1)
    @test t(T(0.5), T(0.5)) == cart(0.5, 0.5)
    @test_throws DomainError((T(-0.5), T(0.0)), "invalid barycentric coordinates for triangle.") t(T(-0.5), T(0.0))
    @test_throws DomainError((T(1), T(1)), "invalid barycentric coordinates for triangle.") t(T(1), T(1))
    @test !hasholes(t)
    @test unique(t) == t
    @test boundary(t) == first(rings(t))
    @test rings(t) == [Ring(cart(0, 0), cart(1, 0), cart(0, 1))]
    @test convexhull(t) == t

    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    equaltest(t)
    isapproxtest(t)

    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test perimeter(t) ≈ T(1 + 1 + √2) * u"m"

    # https://github.com/JuliaGeometry/Meshes.jl/issues/333
    t = Triangle((0.0f0, 0.0f0), (1.0f0, 0.0f0), (0.5f0, 1.0f0))
    @test Point(0.5f0, 0.5f0) ∈ t
    @test Point(0.5e0, 0.5e0) ∈ t

    # circular equality
    t1 = Triangle(T.((1, 1)), T.((2, 2)), T.((3, 3)))
    t2 = Triangle(T.((2, 2)), T.((3, 3)), T.((1, 1)))
    t3 = Triangle(T.((3, 3)), T.((1, 1)), T.((2, 2)))
    @test t1 ≗ t2 ≗ t3

    # point at edge of triangle
    @test cart(3, 1) ∈ Triangle(cart(1, 1), cart(5, 1), cart(3, 3))

    # test angles
    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test all(isapprox.(rad2deg.(angles(t)), T[-90, -45, -45] * u"°", atol=8 * eps(T)))
    @test all(isapprox.(rad2deg.(innerangles(t)), T[90, 45, 45] * u"°", atol=8 * eps(T)))

    # Triangle in 3D space
    t = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0))
    @test area(t) == T(0.5) * u"m^2"
    t = Triangle(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 1))
    @test area(t) > T(0.7) * u"m^2"
    for p in cart.([(0, 0, 0), (1, 0, 0), (0, 1, 1), (0, 0.2, 0.2)])
      @test p ∈ t
    end
    for p in cart.([(-1, 0, 0), (1, 2, 0), (0, 1, 2)])
      @test p ∉ t
    end
    t = Triangle(cart(0, 0, 0), cart(0, 1, 0), cart(0, 0, 1))
    @test t(T(0.0), T(0.0)) == cart(0, 0, 0)
    @test t(T(1.0), T(0.0)) == cart(0, 1, 0)
    @test t(T(0.0), T(1.0)) == cart(0, 0, 1)
    @test t(T(0.5), T(0.5)) == cart(0, 0.5, 0.5)
    @test_throws DomainError((T(-0.5), T(0.0)), "invalid barycentric coordinates for triangle.") t(T(-0.5), T(0.0))
    @test_throws DomainError((T(1), T(1)), "invalid barycentric coordinates for triangle.") t(T(1), T(1))
    @test isapprox(normal(t), vector(1, 0, 0))
    @test isapprox(norm(normal(t)), oneunit(ℳ))
    t = Triangle(cart(0, 0, 0), cart(2, 0, 0), cart(0, 2, 2))
    @test isapprox(normal(t), vector(0, -0.7071067811865475, 0.7071067811865475))
    @test isapprox(norm(normal(t)), oneunit(ℳ))

    # CRS propagation
    t = Triangle(merc(0, 0), merc(1, 0), merc(0, 1))
    @test crs(t(T(0), T(0))) === crs(t)

    t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
    @test sprint(show, t) == "Triangle((x: 0.0 m, y: 0.0 m), (x: 1.0 m, y: 0.0 m), (x: 0.0 m, y: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), t) == """
      Triangle
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m)
      └─ Point(x: 0.0f0 m, y: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), t) == """
      Triangle
      ├─ Point(x: 0.0 m, y: 0.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m)
      └─ Point(x: 0.0 m, y: 1.0 m)"""
    end

    # -----------
    # QUADRANGLE
    # -----------

    # Quadrangle in 2D space
    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test crs(q) <: Cartesian{NoDatum}
    @test Meshes.lentype(q) == ℳ
    @test vertex(q, 1) == cart(0, 0)
    @test vertex(q, 2) == cart(1, 0)
    @test vertex(q, 3) == cart(1, 1)
    @test vertex(q, 4) == cart(0, 1)
    @test area(q) == T(1) * u"m^2"
    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1.5, 1.0), cart(0.5, 1.0))
    @test area(q) == T(1) * u"m^2"
    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1.5, 1.0), cart(0.5, 1.0))
    for p in cart.([(0, 0), (1, 0), (1.5, 1.0), (0.5, 1.0), (0.5, 0.5)])
      @test p ∈ q
    end
    for p in cart.([(0, 1), (1.5, 0.0)])
      @test p ∉ q
    end
    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test !hasholes(q)
    @test unique(q) == q
    @test boundary(q) == first(rings(q))
    @test rings(q) == [Ring(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))]
    @test q(T(0), T(0)) == cart(0, 0)
    @test q(T(1), T(0)) == cart(1, 0)
    @test q(T(1), T(1)) == cart(1, 1)
    @test q(T(0), T(1)) == cart(0, 1)

    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    equaltest(q)
    isapproxtest(q)

    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test_throws DomainError((T(1.2), T(1.2)), "q(u, v) is not defined for u, v outside [0, 1]².") q(T(1.2), T(1.2))

    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test perimeter(q) ≈ T(4) * u"m"

    # Quadrangle in 3D space
    q = Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0))
    @test area(q) == T(1) * u"m^2"
    q = Quadrangle(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 1))
    @test area(q) > T(1) * u"m^2"
    @test q(T(0), T(0)) == cart(0, 0, 0)
    @test q(T(1), T(0)) == cart(1, 0, 0)
    @test q(T(1), T(1)) == cart(1, 1, 0)
    @test q(T(0), T(1)) == cart(0, 1, 1)

    # CRS propagation
    q = Quadrangle(merc(0, 0), merc(1, 0), merc(1, 1), merc(0, 1))
    @test crs(q(T(0), T(0))) === crs(q)

    q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
    @test sprint(show, q) == "Quadrangle((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), q) == """
      Quadrangle
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 1.0f0 m)
      └─ Point(x: 0.0f0 m, y: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), q) == """
      Quadrangle
      ├─ Point(x: 0.0 m, y: 0.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m)
      ├─ Point(x: 1.0 m, y: 1.0 m)
      └─ Point(x: 0.0 m, y: 1.0 m)"""
    end
  end

  @testitem "PolyAreas" begin
    @test paramdim(PolyArea) == 2

    # equality and approximate equality
    outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test poly == poly
    @test poly ≈ poly
    @test crs(poly) <: Cartesian{NoDatum}
    @test Meshes.lentype(poly) == ℳ

    p = PolyArea(cart(0, 0), cart(1, 0), cart(0, 1))
    equaltest(p)
    isapproxtest(p)

    # COMMAND USED TO GENERATE TEST FILES (VARY --seed = 1, 2, ..., 5)
    # rpg --cluster 30 --algo 2opt --format line --seed 1 --output poly1
    fnames = ["poly$i.line" for i in 1:5]
    polys1 = [readpoly(T, joinpath(datadir, fname)) for fname in fnames]
    for poly in polys1
      @test !hasholes(poly)
      @test issimple(poly)
      @test boundary(poly) == first(rings(poly))
      @test nvertices(poly) == 30
      @test orientation(poly) == CCW
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
      @test orientation(poly) == CCW
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
      o = orientation(poly)
      @test o[1] == CCW
      @test all(o[2:end] .== CW)
      @test unique(poly) == poly
    end

    # test bridges
    for poly in [polys1; polys2; polys3]
      b = poly |> Bridge()
      nb = nvertices(b)
      np = nvertices.(rings(poly))
      @test nb ≥ sum(np)
      # orientation always works even
      # in the presence of self-intersections
      @test orientation(b) == CCW
    end

    # test uniqueness
    points = cart.([(1, 1), (2, 2), (2, 2), (3, 3)])
    poly = PolyArea(points)
    unique!(poly)
    @test first(rings(poly)) == Ring(cart.([(1, 1), (2, 2), (3, 3)]))

    # approximately equal vertices
    poly = PolyArea(
      cart.(
        [
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
    )
    upoly = unique(poly)
    @test nvertices(upoly) < nvertices(poly)
    if T === Float32
      @test nvertices(upoly) == 10
    else
      @test nvertices(upoly) == 17
    end

    # invalid inner
    outer = Ring(randpoint2(10))
    p1, p2 = randpoint2(2)
    inner = Ring(p1, p1, p2)
    poly = PolyArea([outer, inner])
    upoly = unique(poly)
    @test hasholes(poly)
    @test !hasholes(upoly)

    # centroid
    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test centroid(poly) == cart(0.5, 0.5)

    # single vertex access
    poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test vertex(poly, 1) == cart(0, 0)
    @test vertex(poly, 4) == cart(0, 1)

    # point in polygonal area
    outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test all(p ∈ poly for p in outer)
    @test cart(0.5, 0.5) ∈ poly
    @test cart(0.2, 0.6) ∈ poly
    @test cart(1.5, 0.5) ∉ poly
    @test cart(-0.5, 0.5) ∉ poly
    @test cart(0.25, 0.25) ∉ poly
    @test cart(0.75, 0.25) ∉ poly
    @test cart(0.75, 0.75) ∈ poly

    # area
    outer = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    hole1 = Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))
    hole2 = Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))
    poly = PolyArea([outer, reverse(hole1), reverse(hole2)])
    @test area(poly) ≈ T(0.92) * u"m^2"

    outer = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    hole1 = Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))
    hole2 = Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))
    poly1 = PolyArea(outer)
    poly2 = PolyArea([outer, reverse(hole1), reverse(hole2)])
    @test sprint(show, poly1) == "PolyArea((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))"
    @test sprint(show, poly2) == "PolyArea(4-Ring, 4-Ring, 4-Ring)"
    @test sprint(show, MIME("text/plain"), poly1) == """
    PolyArea
      outer
      └─ Ring((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))"""
    @test sprint(show, MIME("text/plain"), poly2) == """
    PolyArea
      outer
      └─ Ring((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))
      inner
      ├─ Ring((x: 0.2 m, y: 0.2 m), ..., (x: 0.4 m, y: 0.2 m))
      └─ Ring((x: 0.6 m, y: 0.2 m), ..., (x: 0.8 m, y: 0.2 m))"""

    # should not repeat the first vertex manually
    @test_throws ArgumentError PolyArea(cart.([(0, 0), (0, 0)]))
    @test_throws ArgumentError PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 0)]))
  end

  @testitem "Polyhedra" begin
    @test paramdim(Tetrahedron) == 3
    @test nvertices(Tetrahedron) == 4

    t = Tetrahedron(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1))
    @test crs(t) <: Cartesian{NoDatum}
    @test Meshes.lentype(t) == ℳ
    @test vertex(t, 1) == cart(0, 0, 0)
    @test vertex(t, 2) == cart(1, 0, 0)
    @test vertex(t, 3) == cart(0, 1, 0)
    @test vertex(t, 4) == cart(0, 0, 1)
    @test measure(t) == T(1 / 6) * u"m^3"
    m = boundary(t)
    n = normal.(m)
    @test m isa Mesh
    @test nvertices(m) == 4
    @test nelements(m) == 4
    @test n[1] == vector(0, 0, -1)
    @test n[2] == vector(0, -1, 0)
    @test n[3] == vector(-1, 0, 0)
    @test all(>(T(0) * u"m"), n[4])
    @test t(T(0), T(0), T(0)) ≈ cart(0, 0, 0)
    @test t(T(1), T(0), T(0)) ≈ cart(1, 0, 0)
    @test t(T(0), T(1), T(0)) ≈ cart(0, 1, 0)
    @test t(T(0), T(0), T(1)) ≈ cart(0, 0, 1)
    @test_throws DomainError((T(1), T(1), T(1)), "invalid barycentric coordinates for tetrahedron.") t(T(1), T(1), T(1))

    t = Tetrahedron(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1))
    equaltest(t)
    isapproxtest(t)

    # CRS propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(0), T(0))
    c3 = Cartesian{WGS84Latest}(T(0), T(1), T(0))
    c4 = Cartesian{WGS84Latest}(T(0), T(0), T(1))
    t = Tetrahedron(Point(c1), Point(c2), Point(c3), Point(c4))
    @test crs(t(T(0), T(0), T(0))) === crs(t)

    t = Tetrahedron(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1))
    @test sprint(show, t) == "Tetrahedron((x: 0.0 m, y: 0.0 m, z: 0.0 m), ..., (x: 0.0 m, y: 0.0 m, z: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), t) == """
      Tetrahedron
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 0.0f0 m, y: 1.0f0 m, z: 0.0f0 m)
      └─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), t) == """
      Tetrahedron
      ├─ Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 0.0 m, y: 1.0 m, z: 0.0 m)
      └─ Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)"""
    end

    @test paramdim(Hexahedron) == 3
    @test nvertices(Hexahedron) == 8

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
    @test crs(h) <: Cartesian{NoDatum}
    @test Meshes.lentype(h) == ℳ
    @test vertex(h, 1) == cart(0, 0, 0)
    @test vertex(h, 8) == cart(0, 1, 1)
    @test h(T(0), T(0), T(0)) == cart(0, 0, 0)
    @test h(T(0), T(0), T(1)) == cart(0, 0, 1)
    @test h(T(0), T(1), T(0)) == cart(0, 1, 0)
    @test h(T(0), T(1), T(1)) == cart(0, 1, 1)
    @test h(T(1), T(0), T(0)) == cart(1, 0, 0)
    @test h(T(1), T(0), T(1)) == cart(1, 0, 1)
    @test h(T(1), T(1), T(0)) == cart(1, 1, 0)
    @test h(T(1), T(1), T(1)) == cart(1, 1, 1)

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
    equaltest(h)
    isapproxtest(h)

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
    @test volume(h) ≈ T(1 * 1 * 1) * u"m^3"
    h = Hexahedron(
      cart(0, 0, 0),
      cart(2, 0, 0),
      cart(2, 2, 0),
      cart(0, 2, 0),
      cart(0, 0, 2),
      cart(2, 0, 2),
      cart(2, 2, 2),
      cart(0, 2, 2)
    )
    @test volume(h) ≈ T(2 * 2 * 2) * u"m^3"

    # volume formula of a frustum of a prism is V = 1/3*H*(S₁+S₂+sqrt(S₁*S₂))
    # here we build a hexahedron which is a frustum of a prism with
    # bottom area S₁= 4, top area S₂= 1, height H = 2
    h = Hexahedron(
      cart(0, 0, 0),
      cart(2, 0, 0),
      cart(2, 2, 0),
      cart(0, 2, 0),
      cart(0, 0, 2),
      cart(1, 0, 2),
      cart(1, 1, 2),
      cart(0, 1, 2)
    )
    @test volume(h) ≈ T(1 / 3 * 2 * (1 + 4 + sqrt(1 * 4))) * u"m^3"

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
    m = boundary(h)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    # CRS propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(0), T(0))
    c3 = Cartesian{WGS84Latest}(T(1), T(1), T(0))
    c4 = Cartesian{WGS84Latest}(T(0), T(1), T(0))
    c5 = Cartesian{WGS84Latest}(T(0), T(0), T(1))
    c6 = Cartesian{WGS84Latest}(T(1), T(0), T(1))
    c7 = Cartesian{WGS84Latest}(T(1), T(1), T(1))
    c8 = Cartesian{WGS84Latest}(T(0), T(1), T(1))
    h = Hexahedron(Point(c1), Point(c2), Point(c3), Point(c4), Point(c5), Point(c6), Point(c7), Point(c8))
    @test crs(h(T(0), T(0), T(0))) === crs(h)

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
    @test sprint(show, h) == "Hexahedron((x: 0.0 m, y: 0.0 m, z: 0.0 m), ..., (x: 0.0 m, y: 1.0 m, z: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), h) == """
      Hexahedron
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 1.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 0.0f0 m, y: 1.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 1.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m, z: 1.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 1.0f0 m, z: 1.0f0 m)
      └─ Point(x: 0.0f0 m, y: 1.0f0 m, z: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), h) == """
      Hexahedron
      ├─ Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 1.0 m, y: 1.0 m, z: 0.0 m)
      ├─ Point(x: 0.0 m, y: 1.0 m, z: 0.0 m)
      ├─ Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m, z: 1.0 m)
      ├─ Point(x: 1.0 m, y: 1.0 m, z: 1.0 m)
      └─ Point(x: 0.0 m, y: 1.0 m, z: 1.0 m)"""
    end

    @test paramdim(Pyramid) == 3
    @test nvertices(Pyramid) == 5

    p = Pyramid(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0), cart(0, 0, 1))
    @test crs(p) <: Cartesian{NoDatum}
    @test Meshes.lentype(p) == ℳ
    @test volume(p) ≈ T(1 / 3) * u"m^3"
    m = boundary(p)
    @test m isa Mesh
    @test nelements(m) == 5
    @test m[1] isa Quadrangle
    @test m[2] isa Triangle
    @test m[3] isa Triangle
    @test m[4] isa Triangle
    @test m[5] isa Triangle
    equaltest(p)
    isapproxtest(p)

    p = Pyramid(cart(0, 0, 0), cart(1, 0, 0), cart(1, 1, 0), cart(0, 1, 0), cart(0, 0, 1))
    @test sprint(show, p) == "Pyramid((x: 0.0 m, y: 0.0 m, z: 0.0 m), ..., (x: 0.0 m, y: 0.0 m, z: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), p) == """
      Pyramid
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 1.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 0.0f0 m, y: 1.0f0 m, z: 0.0f0 m)
      └─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), p) == """
      Pyramid
      ├─ Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 1.0 m, y: 1.0 m, z: 0.0 m)
      ├─ Point(x: 0.0 m, y: 1.0 m, z: 0.0 m)
      └─ Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)"""
    end

    @test paramdim(Wedge) == 3
    @test nvertices(Wedge) == 6

    w = Wedge(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1), cart(1, 0, 1), cart(0, 1, 1))
    @test crs(w) <: Cartesian{NoDatum}
    @test Meshes.lentype(w) == ℳ
    @test volume(w) ≈ T(1 / 2) * u"m^3"
    m = boundary(w)
    @test m isa Mesh
    @test nelements(m) == 5
    @test m[1] isa Triangle
    @test m[2] isa Triangle
    @test m[3] isa Quadrangle
    @test m[4] isa Quadrangle
    @test m[5] isa Quadrangle
    equaltest(w)
    isapproxtest(w)

    w = Wedge(cart(0, 0, 0), cart(1, 0, 0), cart(0, 1, 0), cart(0, 0, 1), cart(1, 0, 1), cart(0, 1, 1))
    @test sprint(show, w) == "Wedge((x: 0.0 m, y: 0.0 m, z: 0.0 m), ..., (x: 0.0 m, y: 1.0 m, z: 1.0 m))"
    if T === Float32
      @test sprint(show, MIME("text/plain"), w) == """
      Wedge
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 0.0f0 m, y: 1.0f0 m, z: 0.0f0 m)
      ├─ Point(x: 0.0f0 m, y: 0.0f0 m, z: 1.0f0 m)
      ├─ Point(x: 1.0f0 m, y: 0.0f0 m, z: 1.0f0 m)
      └─ Point(x: 0.0f0 m, y: 1.0f0 m, z: 1.0f0 m)"""
    else
      @test sprint(show, MIME("text/plain"), w) == """
      Wedge
      ├─ Point(x: 0.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m, z: 0.0 m)
      ├─ Point(x: 0.0 m, y: 1.0 m, z: 0.0 m)
      ├─ Point(x: 0.0 m, y: 0.0 m, z: 1.0 m)
      ├─ Point(x: 1.0 m, y: 0.0 m, z: 1.0 m)
      └─ Point(x: 0.0 m, y: 1.0 m, z: 1.0 m)"""
    end
  end
end
