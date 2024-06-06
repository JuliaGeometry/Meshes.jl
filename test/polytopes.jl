@testset "Polytopes" begin
  @testset "Segments" begin
    @test paramdim(Segment) == 1
    @test nvertices(Segment) == 2

    s = Segment(point(1.0), point(2.0))
    @test Meshes.crs(s) <: Cartesian{NoDatum}
    @test Meshes.lentype(s) == ℳ
    @test vertex(s, 1) == point(1.0)
    @test vertex(s, 2) == point(2.0)
    @test all(point(x) ∈ s for x in 1:0.01:2)
    @test all(point(x) ∉ s for x in [-1.0, 0.0, 0.99, 2.1, 5.0, 10.0])
    @test s ≈ s
    @test !(s ≈ Segment(point(2.0), point(1.0)))
    @test !(s ≈ Segment(point(-1.0), point(2.0)))

    s = Segment(point(0, 0), point(1, 1))
    @test minimum(s) == point(0, 0)
    @test maximum(s) == point(1, 1)
    @test extrema(s) == (point(0, 0), point(1, 1))
    @test isapprox(length(s), sqrt(T(2)) * u"m")
    @test s(T(0)) == point(0, 0)
    @test s(T(1)) == point(1, 1)
    @test all(point(x, x) ∈ s for x in 0:0.01:1)
    @test all(p ∉ s for p in [point(-0.1, -0.1), point(1.1, 1.1), point(0.5, 0.49), point(1, 2)])
    @test_throws DomainError(T(1.2), "s(t) is not defined for t outside [0, 1].") s(T(1.2))
    @test_throws DomainError(T(-0.5), "s(t) is not defined for t outside [0, 1].") s(T(-0.5))
    @test s ≈ s
    @test !(s ≈ Segment(point(1, 1), point(0, 0)))
    @test !(s ≈ Segment(point(1, 2), point(0, 0)))

    s = Segment(point(0, 0, 0), point(1, 1, 1))
    @test all(point(x, x, x) ∈ s for x in 0:0.01:1)
    @test all(p ∉ s for p in [point(-0.1, -0.1, -0.1), point(1.1, 1.1, 1.1)])
    @test all(p ∉ s for p in [point(0.5, 0.5, 0.49), point(1, 1, 2)])
    @test s ≈ s
    @test !(s ≈ Segment(point(1, 1, 1), point(0, 0, 0)))
    @test !(s ≈ Segment(point(1, 1, 1), point(0, 1, 0)))

    s = Segment(Point(1.0, 1.0, 1.0, 1.0), Point(2.0, 2.0, 2.0, 2.0))
    @test all(Point(x, x, x, x) ∈ s for x in 1:0.01:2)
    @test all(p ∉ s for p in [Point(0.99, 0.99, 0.99, 0.99), Point(2.1, 2.1, 2.1, 2.1)])
    @test all(p ∉ s for p in [Point(1.5, 1.5, 1.5, 1.49), Point(1, 1, 2, 1.0)])
    @test s ≈ s
    @test !(s ≈ Segment(Point(2, 2, 2, 2), Point(1, 1, 1, 1)))
    @test !(s ≈ Segment(Point(1, 1, 2, 1), Point(0, 0, 0, 0)))

    s = Segment(point(0, 0, 0), point(1, 1, 1))
    @test boundary(s) == Multi([point(0, 0, 0), point(1, 1, 1)])
    @test perimeter(s) == zero(T) * u"m"
    @test center(s) == point(0.5, 0.5, 0.5)
    @test Meshes.lentype(center(s)) == ℳ

    # unitful coordinates
    x1 = T(0)u"m"
    x2 = T(1)u"m"
    s = Segment(Point(x1, x1, x1), Point(x2, x2, x2))
    @test boundary(s) == Multi([Point(x1, x1, x1), Point(x2, x2, x2)])
    @test perimeter(s) == 0u"m"
    xm = T(0.5)u"m"
    @test center(s) == Point(xm, xm, xm)
    @test Meshes.lentype(center(s)) == typeof(xm)

    s = rand(Segment{2})
    @test s isa Segment
    @test embeddim(s) == 2
    @test Meshes.lentype(s) === Meshes.Met{Float64}
    s = rand(Segment{3})
    @test s isa Segment
    @test embeddim(s) == 3
    @test Meshes.lentype(s) === Meshes.Met{Float64}

    # datum propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(1))
    s = Segment(Point(c1), Point(c2))
    @test datum(Meshes.crs(s(T(0)))) === WGS84Latest

    s = Segment(point(0, 0), point(1, 1))
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

  @testset "Ropes/Rings" begin
    c1 = Rope(point.([(1, 1), (2, 2)]))
    c2 = Rope(point(1, 1), point(2, 2))
    c3 = Rope(T.((1, 1.0)), T.((2.0, 2.0)))
    @test c1 == c2 == c3
    c1 = Ring(point.([(1, 1), (2, 2)]))
    c2 = Ring(point(1, 1), point(2, 2))
    c3 = Ring(T.((1, 1.0)), T.((2.0, 2.0)))
    @test c1 == c2 == c3

    c = Rope(point.([(1, 1), (2, 2)]))
    @test Meshes.crs(c) <: Cartesian{NoDatum}
    @test Meshes.lentype(c) == ℳ
    @test vertex(c, 1) == point(1, 1)
    @test vertex(c, 2) == point(2, 2)
    c = Ring(point.([(1, 1), (2, 2)]))
    @test Meshes.crs(c) <: Cartesian{NoDatum}
    @test Meshes.lentype(c) == ℳ
    @test vertex(c, 0) == point(2, 2)
    @test vertex(c, 1) == point(1, 1)
    @test vertex(c, 2) == point(2, 2)
    @test vertex(c, 3) == point(1, 1)
    @test vertex(c, 4) == point(2, 2)

    c = Rope(point.([(1, 1), (2, 2), (3, 3)]))
    @test collect(segments(c)) == [Segment(point(1, 1), point(2, 2)), Segment(point(2, 2), point(3, 3))]
    c = Ring(point.([(1, 1), (2, 2), (3, 3)]))
    @test collect(segments(c)) ==
          [Segment(point(1, 1), point(2, 2)), Segment(point(2, 2), point(3, 3)), Segment(point(3, 3), point(1, 1))]

    c = Rope(point.([(1, 1), (2, 2), (2, 2), (3, 3)]))
    @test unique(c) == Rope(point.([(1, 1), (2, 2), (3, 3)]))
    @test c == Rope(point.([(1, 1), (2, 2), (2, 2), (3, 3)]))
    unique!(c)
    @test c == Rope(point.([(1, 1), (2, 2), (3, 3)]))

    c = Rope(point.([(1, 1), (2, 2), (3, 3)]))
    @test close(c) == Ring(point.([(1, 1), (2, 2), (3, 3)]))
    c = Ring(point.([(1, 1), (2, 2), (3, 3)]))
    @test open(c) == Rope(point.([(1, 1), (2, 2), (3, 3)]))

    c = Rope(point.([(1, 1), (2, 2), (3, 3)]))
    reverse!(c)
    @test c == Rope(point.([(3, 3), (2, 2), (1, 1)]))
    c = Rope(point.([(1, 1), (2, 2), (3, 3)]))
    @test reverse(c) == Rope(point.([(3, 3), (2, 2), (1, 1)]))

    c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test angles(c) ≈ [-π / 2, -π / 2, -π / 2, -π / 2]
    c = Rope(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test angles(c) ≈ [-π / 2, -π / 2]
    c = Ring(point.([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2)]))
    @test angles(c) ≈ [-atan(2), -π / 2, +π / 2, -π / 2, -π / 2, -(π - atan(2))]
    @test innerangles(c) ≈ [atan(2), π / 2, 3π / 2, π / 2, π / 2, π - atan(2)]

    c1 = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    c2 = Ring(vertices(c1))
    @test c1 == c2

    c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test centroid(c) == point(0.5, 0.5)

    c = Rope(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test boundary(c) == Multi(point.([(0, 0), (0, 1)]))
    c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test isnothing(boundary(c))

    # should not repeat the first vertex manually
    @test_throws ArgumentError Ring(point.([(0, 0), (0, 0)]))
    @test_throws ArgumentError Ring(point.([(0, 0), (1, 0), (1, 1), (0, 0)]))

    # degenerate rings with 1 or 2 vertices are allowed
    r = Ring(point.([(0, 0)]))
    @test isclosed(r)
    @test nvertices(r) == 1
    @test collect(segments(r)) == [Segment(point(0, 0), point(0, 0))]
    r = Ring(point.([(0, 0), (1, 1)]))
    @test isclosed(r)
    @test nvertices(r) == 2
    @test collect(segments(r)) == [Segment(point(0, 0), point(1, 1)), Segment(point(1, 1), point(0, 0))]

    p1 = point(1, 1)
    p2 = point(3, 1)
    p3 = point(1, 0)
    p4 = point(3, 0)
    pts = point.([(0, 0), (2, 2), (4, 0)])
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
      point.(
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

    r = rand(Rope{2})
    @test r isa Rope
    @test embeddim(r) == 2
    @test Meshes.lentype(r) === Meshes.Met{Float64}
    r = rand(Rope{3})
    @test r isa Rope
    @test embeddim(r) == 3
    @test Meshes.lentype(r) === Meshes.Met{Float64}

    r = rand(Ring{2})
    @test r isa Ring
    @test embeddim(r) == 2
    @test Meshes.lentype(r) === Meshes.Met{Float64}
    r = rand(Ring{3})
    @test r isa Ring
    @test embeddim(r) == 3
    @test Meshes.lentype(r) === Meshes.Met{Float64}

    # issimple benchmark
    r = Sphere(point(0, 0), T(1)) |> pointify |> Ring
    @test issimple(r)
    @test @elapsed(issimple(r)) < 0.02
    @test @allocated(issimple(r)) < 950000

    # innerangles in 3D is obtained via projection
    r1 = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    r2 = Ring(point.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
    @test innerangles(r1) ≈ innerangles(r2)

    # datum propagation
    tuples = [T.((0, 0)), T.((1, 0)), T.((1, 1)), T.((0, 1))]
    points = Point.(Cartesian{WGS84Latest}.(tuples))
    r = Ring(points)
    @test datum(Meshes.crs(centroid(r))) === WGS84Latest

    ri = Ring(point.([(1, 1), (2, 2), (3, 3)]))
    ro = Rope(point.([(1, 1), (2, 2), (3, 3)]))
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

  @testset "Ngons" begin
    pts = (point(0, 0), point(1, 0), point(0, 1))
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

      n = rand(NGON{2})
      @test n isa NGON
      @test embeddim(n) == 2
      @test Meshes.lentype(n) === Meshes.Met{Float64}
      n = rand(NGON{3})
      @test n isa NGON
      @test embeddim(n) == 3
      @test Meshes.lentype(n) === Meshes.Met{Float64}
    end

    # error: the number of vertices must be greater than or equal to 3
    @test_throws ArgumentError Ngon(point(0, 0), point(1, 1))
    @test_throws ArgumentError Ngon{2}(point(0, 0), point(1, 1))

    # ---------
    # TRIANGLE
    # ---------

    # Triangle in 2D space
    t = Triangle(point(0, 0), point(1, 0), point(0, 1))
    @test Meshes.crs(t) <: Cartesian{NoDatum}
    @test Meshes.lentype(t) == ℳ
    @test vertex(t, 1) == point(0, 0)
    @test vertex(t, 2) == point(1, 0)
    @test vertex(t, 3) == point(0, 1)
    @test signarea(t) == T(0.5) * u"m^2"
    @test area(t) == T(0.5) * u"m^2"
    t = Triangle(point(0, 0), point(0, 1), point(1, 0))
    @test signarea(t) == T(-0.5) * u"m^2"
    @test area(t) == T(0.5) * u"m^2"
    t = Triangle(point(0, 0), point(1, 0), point(1, 1))
    for p in point.([(0, 0), (1, 0), (1, 1), (0.5, 0.0), (1.0, 0.5), (0.5, 0.5)])
      @test p ∈ t
    end
    for p in point.([(-1, 0), (0, -1), (0.5, 1.0)])
      @test p ∉ t
    end
    t = Triangle(point(0.4, 0.4), point(0.6, 0.4), point(0.8, 0.4))
    @test point(0.2, 0.4) ∉ t
    t = Triangle(point(0, 0), point(1, 0), point(0, 1))
    @test t(T(0.0), T(0.0)) == point(0, 0)
    @test t(T(1.0), T(0.0)) == point(1, 0)
    @test t(T(0.0), T(1.0)) == point(0, 1)
    @test t(T(0.5), T(0.5)) == point(0.5, 0.5)
    @test_throws DomainError((T(-0.5), T(0.0)), "invalid barycentric coordinates for triangle.") t(T(-0.5), T(0.0))
    @test_throws DomainError((T(1), T(1)), "invalid barycentric coordinates for triangle.") t(T(1), T(1))
    @test !hasholes(t)
    @test unique(t) == t
    @test boundary(t) == first(rings(t))
    @test rings(t) == [Ring(point(0, 0), point(1, 0), point(0, 1))]
    @test convexhull(t) == t

    t = Triangle(point(0, 0), point(1, 0), point(0, 1))
    @test perimeter(t) ≈ T(1 + 1 + √2) * u"m"

    # https://github.com/JuliaGeometry/Meshes.jl/issues/333
    t = Triangle((0.0f0, 0.0f0), (1.0f0, 0.0f0), (0.5f0, 1.0f0))
    @test Point(0.5f0, 0.5f0) ∈ t
    @test Point(0.5e0, 0.5e0) ∈ t

    # point at edge of triangle
    @test point(3, 1) ∈ Triangle(point(1, 1), point(5, 1), point(3, 3))

    # test angles
    t = Triangle(point(0, 0), point(1, 0), point(0, 1))
    @test all(isapprox.(rad2deg.(angles(t)), T[-90, -45, -45], atol=8 * eps(T)))
    @test all(isapprox.(rad2deg.(innerangles(t)), T[90, 45, 45], atol=8 * eps(T)))

    # Triangle in 3D space
    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))
    @test area(t) == T(0.5) * u"m^2"
    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 1))
    @test area(t) > T(0.7) * u"m^2"
    for p in point.([(0, 0, 0), (1, 0, 0), (0, 1, 1), (0, 0.2, 0.2)])
      @test p ∈ t
    end
    for p in point.([(-1, 0, 0), (1, 2, 0), (0, 1, 2)])
      @test p ∉ t
    end
    t = Triangle(point(0, 0, 0), point(0, 1, 0), point(0, 0, 1))
    @test t(T(0.0), T(0.0)) == point(0, 0, 0)
    @test t(T(1.0), T(0.0)) == point(0, 1, 0)
    @test t(T(0.0), T(1.0)) == point(0, 0, 1)
    @test t(T(0.5), T(0.5)) == point(0, 0.5, 0.5)
    @test_throws DomainError((T(-0.5), T(0.0)), "invalid barycentric coordinates for triangle.") t(T(-0.5), T(0.0))
    @test_throws DomainError((T(1), T(1)), "invalid barycentric coordinates for triangle.") t(T(1), T(1))
    @test isapprox(normal(t), vector(1, 0, 0))
    @test isapprox(norm(normal(t)), oneunit(ℳ))
    t = Triangle(point(0, 0, 0), point(2, 0, 0), point(0, 2, 2))
    @test isapprox(normal(t), vector(0, -0.7071067811865475, 0.7071067811865475))
    @test isapprox(norm(normal(t)), oneunit(ℳ))

    t = Triangle(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0))
    @test_throws ErrorException("signed area only defined for triangles embedded in R², use `area` instead") signarea(t)

    # datum propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(0))
    c3 = Cartesian{WGS84Latest}(T(0), T(1))
    t = Triangle(Point(c1), Point(c2), Point(c3))
    @test datum(Meshes.crs(t(T(0), T(0)))) === WGS84Latest

    t = Triangle(point(0, 0), point(1, 0), point(0, 1))
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

    # test periodicity of Quadrangle
    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))

    # Quadrangle in 2D space
    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    @test Meshes.crs(q) <: Cartesian{NoDatum}
    @test Meshes.lentype(q) == ℳ
    @test vertex(q, 1) == point(0, 0)
    @test vertex(q, 2) == point(1, 0)
    @test vertex(q, 3) == point(1, 1)
    @test vertex(q, 4) == point(0, 1)
    @test area(q) == T(1) * u"m^2"
    q = Quadrangle(point(0, 0), point(1, 0), point(1.5, 1.0), point(0.5, 1.0))
    @test area(q) == T(1) * u"m^2"
    q = Quadrangle(point(0, 0), point(1, 0), point(1.5, 1.0), point(0.5, 1.0))
    for p in point.([(0, 0), (1, 0), (1.5, 1.0), (0.5, 1.0), (0.5, 0.5)])
      @test p ∈ q
    end
    for p in point.([(0, 1), (1.5, 0.0)])
      @test p ∉ q
    end
    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    @test !hasholes(q)
    @test unique(q) == q
    @test boundary(q) == first(rings(q))
    @test rings(q) == [Ring(point(0, 0), point(1, 0), point(1, 1), point(0, 1))]
    @test q(T(0), T(0)) == point(0, 0)
    @test q(T(1), T(0)) == point(1, 0)
    @test q(T(1), T(1)) == point(1, 1)
    @test q(T(0), T(1)) == point(0, 1)

    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    @test_throws DomainError((T(1.2), T(1.2)), "q(u, v) is not defined for u, v outside [0, 1]².") q(T(1.2), T(1.2))

    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
    @test perimeter(q) ≈ T(4) * u"m"

    # Quadrangle in 3D space
    q = Quadrangle(point(0, 0, 0), point(1, 0, 0), point(1, 1, 0), point(0, 1, 0))
    @test area(q) == T(1) * u"m^2"
    q = Quadrangle(point(0, 0, 0), point(1, 0, 0), point(1, 1, 0), point(0, 1, 1))
    @test area(q) > T(1) * u"m^2"
    @test q(T(0), T(0)) == point(0, 0, 0)
    @test q(T(1), T(0)) == point(1, 0, 0)
    @test q(T(1), T(1)) == point(1, 1, 0)
    @test q(T(0), T(1)) == point(0, 1, 1)

    # datum propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(0))
    c3 = Cartesian{WGS84Latest}(T(1), T(1))
    c4 = Cartesian{WGS84Latest}(T(0), T(1))
    q = Quadrangle(Point(c1), Point(c2), Point(c3), Point(c4))
    @test datum(Meshes.crs(q(T(0), T(0)))) === WGS84Latest

    q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
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

  @testset "PolyAreas" begin
    @test paramdim(PolyArea) == 2

    # equality and approximate equality
    outer = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test poly == poly
    @test poly ≈ poly
    @test Meshes.crs(poly) <: Cartesian{NoDatum}
    @test Meshes.lentype(poly) == ℳ

    # outer chain with 2 vertices is fixed by default
    poly = PolyArea(point.([(0, 0), (1, 0)]))
    @test rings(poly) == [Ring(point.([(0, 0), (0.5, 0.0), (1, 0)]))]

    # inner chain with 2 vertices is removed by default
    poly = PolyArea([point.([(0, 0), (1, 0), (1, 1), (0, 1)]), point.([(1, 2), (2, 3)])])
    @test rings(poly) == [Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))]

    # orientation of chains is fixed by default
    poly = PolyArea(point.([(0, 0), (0, 1), (1, 1), (1, 0)]))
    @test vertices(poly) == CircularVector(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    poly = PolyArea(point.([(0, 0), (0, 1), (1, 1), (1, 0)]), fix=false)
    @test vertices(poly) == CircularVector(point.([(0, 0), (0, 1), (1, 1), (1, 0)]))

    # test accessor methods
    poly = PolyArea(point.([(1, 2), (2, 3)]), fix=false)
    @test vertices(poly) == CircularVector(point.([(1, 2), (2, 3)]))
    poly = PolyArea([point.([(1, 2), (2, 3)]), point.([(1.1, 2.54), (1.4, 1.5)])], fix=false)
    @test vertices(poly) == CircularVector(point.([(1, 2), (2, 3), (1.1, 2.54), (1.4, 1.5)]))

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
        @test orientation(poly, algo) == CCW
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
        @test orientation(poly, algo) == CCW
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
        @test orients[1] == CCW
        @test all(orients[2:end] .== CW)
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
      @test orientation(b, TriangleOrientation()) == CCW
      # winding orientation is only suitable
      # for simple polygonal chains
      # if issimple(b)
      #   @test orientation(b, WindingOrientation()) == CCW
      # end
    end

    # test uniqueness
    points = point.([(1, 1), (2, 2), (2, 2), (3, 3)])
    poly = PolyArea(points)
    unique!(poly)
    @test first(rings(poly)) == Ring(point.([(1, 1), (2, 2), (3, 3)]))

    # approximately equal vertices
    poly = PolyArea(
      point.(
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
    poly = PolyArea(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test centroid(poly) == point(0.5, 0.5)

    # single vertex access
    poly = PolyArea(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
    @test vertex(poly, 1) == point(0, 0)
    @test vertex(poly, 4) == point(0, 1)

    # point in polygonal area
    outer = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test all(p ∈ poly for p in outer)
    @test point(0.5, 0.5) ∈ poly
    @test point(0.2, 0.6) ∈ poly
    @test point(1.5, 0.5) ∉ poly
    @test point(-0.5, 0.5) ∉ poly
    @test point(0.25, 0.25) ∉ poly
    @test point(0.75, 0.25) ∉ poly
    @test point(0.75, 0.75) ∈ poly

    # area
    outer = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly = PolyArea([outer, hole1, hole2])
    @test area(poly) ≈ T(0.92) * u"m^2"

    p = rand(PolyArea{2})
    @test p isa PolyArea
    @test embeddim(p) == 2
    @test Meshes.lentype(p) === Meshes.Met{Float64}
    p = rand(PolyArea{3})
    @test p isa PolyArea
    @test embeddim(p) == 3
    @test Meshes.lentype(p) === Meshes.Met{Float64}

    outer = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
    hole1 = point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
    hole2 = point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
    poly1 = PolyArea(outer)
    poly2 = PolyArea([outer, hole1, hole2])
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
    @test_throws ArgumentError PolyArea(point.([(0, 0), (0, 0)]))
    @test_throws ArgumentError PolyArea(point.([(0, 0), (1, 0), (1, 1), (0, 0)]))
  end

  @testset "Polyhedra" begin
    @test paramdim(Tetrahedron) == 3
    @test nvertices(Tetrahedron) == 4

    t = Tetrahedron(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0), point(0, 0, 1))
    @test Meshes.crs(t) <: Cartesian{NoDatum}
    @test Meshes.lentype(t) == ℳ
    @test vertex(t, 1) == point(0, 0, 0)
    @test vertex(t, 2) == point(1, 0, 0)
    @test vertex(t, 3) == point(0, 1, 0)
    @test vertex(t, 4) == point(0, 0, 1)
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
    @test t(T(0), T(0), T(0)) ≈ point(0, 0, 0)
    @test t(T(1), T(0), T(0)) ≈ point(1, 0, 0)
    @test t(T(0), T(1), T(0)) ≈ point(0, 1, 0)
    @test t(T(0), T(0), T(1)) ≈ point(0, 0, 1)
    @test_throws DomainError((T(1), T(1), T(1)), "invalid barycentric coordinates for tetrahedron.") t(T(1), T(1), T(1))

    t = rand(Tetrahedron{3})
    @test t isa Tetrahedron
    @test embeddim(t) == 3
    @test Meshes.lentype(t) === Meshes.Met{Float64}

    # datum propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(0), T(0))
    c3 = Cartesian{WGS84Latest}(T(0), T(1), T(0))
    c4 = Cartesian{WGS84Latest}(T(0), T(0), T(1))
    t = Tetrahedron(Point(c1), Point(c2), Point(c3), Point(c4))
    @test datum(Meshes.crs(t(T(0), T(0), T(0)))) === WGS84Latest

    t = Tetrahedron(point(0, 0, 0), point(1, 0, 0), point(0, 1, 0), point(0, 0, 1))
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
      point(0, 0, 0),
      point(1, 0, 0),
      point(1, 1, 0),
      point(0, 1, 0),
      point(0, 0, 1),
      point(1, 0, 1),
      point(1, 1, 1),
      point(0, 1, 1)
    )
    @test Meshes.crs(h) <: Cartesian{NoDatum}
    @test Meshes.lentype(h) == ℳ
    @test vertex(h, 1) == point(0, 0, 0)
    @test vertex(h, 8) == point(0, 1, 1)
    @test h(T(0), T(0), T(0)) == point(0, 0, 0)
    @test h(T(0), T(0), T(1)) == point(0, 0, 1)
    @test h(T(0), T(1), T(0)) == point(0, 1, 0)
    @test h(T(0), T(1), T(1)) == point(0, 1, 1)
    @test h(T(1), T(0), T(0)) == point(1, 0, 0)
    @test h(T(1), T(0), T(1)) == point(1, 0, 1)
    @test h(T(1), T(1), T(0)) == point(1, 1, 0)
    @test h(T(1), T(1), T(1)) == point(1, 1, 1)

    h = Hexahedron(
      point(0, 0, 0),
      point(1, 0, 0),
      point(1, 1, 0),
      point(0, 1, 0),
      point(0, 0, 1),
      point(1, 0, 1),
      point(1, 1, 1),
      point(0, 1, 1)
    )
    @test volume(h) ≈ T(1 * 1 * 1) * u"m^3"
    h = Hexahedron(
      point(0, 0, 0),
      point(2, 0, 0),
      point(2, 2, 0),
      point(0, 2, 0),
      point(0, 0, 2),
      point(2, 0, 2),
      point(2, 2, 2),
      point(0, 2, 2)
    )
    @test volume(h) ≈ T(2 * 2 * 2) * u"m^3"

    # volume formula of a frustum of a prism is V = 1/3*H*(S₁+S₂+sqrt(S₁*S₂))
    # here we build a hexahedron which is a frustum of a prism with
    # bottom area S₁= 4, top area S₂= 1, height H = 2
    h = Hexahedron(
      point(0, 0, 0),
      point(2, 0, 0),
      point(2, 2, 0),
      point(0, 2, 0),
      point(0, 0, 2),
      point(1, 0, 2),
      point(1, 1, 2),
      point(0, 1, 2)
    )
    @test volume(h) ≈ T(1 / 3 * 2 * (1 + 4 + sqrt(1 * 4))) * u"m^3"

    h = Hexahedron(
      point(0, 0, 0),
      point(1, 0, 0),
      point(1, 1, 0),
      point(0, 1, 0),
      point(0, 0, 1),
      point(1, 0, 1),
      point(1, 1, 1),
      point(0, 1, 1)
    )
    m = boundary(h)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    h = rand(Hexahedron{3})
    @test h isa Hexahedron
    @test embeddim(h) == 3
    @test Meshes.lentype(h) === Meshes.Met{Float64}

    # datum propagation
    c1 = Cartesian{WGS84Latest}(T(0), T(0), T(0))
    c2 = Cartesian{WGS84Latest}(T(1), T(0), T(0))
    c3 = Cartesian{WGS84Latest}(T(1), T(1), T(0))
    c4 = Cartesian{WGS84Latest}(T(0), T(1), T(0))
    c5 = Cartesian{WGS84Latest}(T(0), T(0), T(1))
    c6 = Cartesian{WGS84Latest}(T(1), T(0), T(1))
    c7 = Cartesian{WGS84Latest}(T(1), T(1), T(1))
    c8 = Cartesian{WGS84Latest}(T(0), T(1), T(1))
    h = Hexahedron(Point(c1), Point(c2), Point(c3), Point(c4), Point(c5), Point(c6), Point(c7), Point(c8))
    @test datum(Meshes.crs(h(T(0), T(0), T(0)))) === WGS84Latest

    h = Hexahedron(
      point(0, 0, 0),
      point(1, 0, 0),
      point(1, 1, 0),
      point(0, 1, 0),
      point(0, 0, 1),
      point(1, 0, 1),
      point(1, 1, 1),
      point(0, 1, 1)
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

    p = Pyramid(point(0, 0, 0), point(1, 0, 0), point(1, 1, 0), point(0, 1, 0), point(0, 0, 1))
    @test Meshes.crs(p) <: Cartesian{NoDatum}
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

    p = rand(Pyramid{3})
    @test p isa Pyramid
    @test embeddim(p) == 3
    @test Meshes.lentype(p) === Meshes.Met{Float64}

    p = Pyramid(point(0, 0, 0), point(1, 0, 0), point(1, 1, 0), point(0, 1, 0), point(0, 0, 1))
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
  end
end
