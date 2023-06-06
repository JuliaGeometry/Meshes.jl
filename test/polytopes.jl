@testset "Polytopes" begin
  @testset "Segments" begin
    @test paramdim(Segment) == 1
    @test nvertices(Segment) == 2
    @test isconvex(Segment)
    @test isperiodic(Segment) == (false,)

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
    @test center(s) == Point(0.5, 0.5, 0.5)
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
  end

  @testset "Ngons" begin
    @test paramdim(Ngon) == 2
    NGONS = [Triangle, Quadrangle, Pentagon, Hexagon, Heptagon, Octagon, Nonagon, Decagon]
    NVERT = 3:10
    for i in 1:length(NGONS)
      @test paramdim(NGONS[i]) == 2
      @test nvertices(NGONS[i]) == NVERT[i]
    end

    # Triangle in 2D space
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
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
    @test bridge(t) == (first(rings(t)), [])

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
    t = Triangle(P2[(0, 0), (1, 0), (0, 1)])
    @test isconvex(t)
    t = Triangle(P3[(0, 0, 0), (1, 0, 0), (0, 1, 0)])
    @test isconvex(t)

    # test periodicity of Quadrangle
    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    @test isperiodic(Quadrangle) == (false, false)
    @test isperiodic(q) == (false, false)

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
    @test bridge(q) == (first(rings(q)), [])
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
  end

  @testset "PolyAreas" begin
    @test paramdim(PolyArea) == 2

    # degenerate outer chain with 2 vertices is allowed
    poly = PolyArea(P2[(0, 0), (1, 0)])
    @test rings(poly) == [Ring(P2[(0, 0), (1, 0)])]

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
      b, _ = bridge(poly)
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
    chain, _ = bridge(poly)
    target = T[
      0.0 0.2 0.2 0.4 0.4 0.6 0.6 0.8 0.8 0.6 0.4 0.2 0.0 1.0 1.0 0.0
      0.0 0.2 0.4 0.4 0.2 0.2 0.4 0.4 0.2 0.2 0.2 0.2 0.0 0.0 1.0 1.0
    ]
    @test vertices(chain) == Point.(Tuple.(eachcol(target)))

    # test uniqueness
    points = P2[(1, 1), (2, 2), (2, 2), (3, 3)]
    poly = PolyArea(points)
    unique!(poly)
    @test first(rings(poly)) == Ring(P2[(1, 1), (2, 2), (3, 3)])

    # unique and bridges
    poly = PolyArea(P2[(0, 0), (1, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1), (0, 1)])
    chain, _ = poly |> unique |> bridge
    @test chain == Ring(P2[(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)])

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
    point = P2(0.5, 0.5)
    bytes = @allocated point ∈ poly
    @test bytes == 0

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

    # should not repeat the first vertex manually
    @test_throws ArgumentError PolyArea(P2[(0, 0), (0, 0)])
    @test_throws ArgumentError PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 0)])
  end

  @testset "Polyhedra" begin
    @test paramdim(Tetrahedron) == 3
    @test nvertices(Tetrahedron) == 4

    t = Tetrahedron(P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
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

    @test paramdim(Hexahedron) == 3
    @test nvertices(Hexahedron) == 8
    @test isperiodic(Hexahedron) == (false, false, false)

    h = Hexahedron(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    @test vertex(h, 1) == P3(0, 0, 0)
    @test vertex(h, 8) == P3(0, 1, 1)
    @test isperiodic(h) == (false, false, false)
    @test h(T(0), T(0), T(0)) == P3(0, 0, 0)
    @test h(T(0), T(0), T(1)) == P3(0, 0, 1)
    @test h(T(0), T(1), T(0)) == P3(0, 1, 0)
    @test h(T(0), T(1), T(1)) == P3(0, 1, 1)
    @test h(T(1), T(0), T(0)) == P3(1, 0, 0)
    @test h(T(1), T(0), T(1)) == P3(1, 0, 1)
    @test h(T(1), T(1), T(0)) == P3(1, 1, 0)
    @test h(T(1), T(1), T(1)) == P3(1, 1, 1)

    h = Hexahedron(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    @test volume(h) ≈ T(1 * 1 * 1)
    h = Hexahedron(P3[(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0), (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2)])
    @test volume(h) ≈ T(2 * 2 * 2)

    # volume formula of a frustum of a prism is V = 1/3*H*(S₁+S₂+sqrt(S₁*S₂))
    # here we build a hexahedron which is a frustum of a prism with
    # bottom area S₁= 4, top area S₂= 1, height H = 2
    h = Hexahedron(P3[(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0), (0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2)])
    @test volume(h) ≈ T(1 / 3 * 2 * (1 + 4 + sqrt(1 * 4)))

    h = Hexahedron(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    m = boundary(h)
    @test m isa Mesh
    @test nvertices(m) == 8
    @test nelements(m) == 6

    @test paramdim(Pyramid) == 3
    @test nvertices(Pyramid) == 5

    p = Pyramid(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1)])
    @test volume(p) ≈ T(1 / 3)
    m = boundary(p)
    @test m isa Mesh
    @test nelements(m) == 5
    @test m[1] isa Quadrangle
    @test m[2] isa Triangle
    @test m[3] isa Triangle
    @test m[4] isa Triangle
    @test m[5] isa Triangle
  end
end
