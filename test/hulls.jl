@testitem "Hulls" setup = [Setup] begin
  for method in [GrahamScan(), JarvisMarch()]
    # basic test
    pts = [cart(rand(T), rand(T)) for _ in 1:10]
    chul = hull(pts, method)
    @test all(pts .∈ Ref(chul))

    # duplicated points
    pts = [cart(rand(T), rand(T)) for _ in 1:10]
    dup = [pts; pts]
    chul = hull(dup, method)
    @test all(pts .∈ Ref(chul))

    # corner cases
    pts = cart.([(0, 0)])
    chul = hull(pts, method)
    @test chul == cart(0, 0)
    pts = cart.([(0, 1), (1, 0)])
    chul = hull(pts, method)
    @test chul == Segment(cart(0, 1), cart(1, 0))
    pts = cart.([(1, 0), (0, 0), (0, 1)])
    chul = hull(pts, method)
    @test vertices(chul) == cart.([(0, 0), (1, 0), (0, 1)])

    # original point set is already in hull
    pts = cart.([(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)])
    chul = hull(pts, method)
    verts = vertices(chul)
    @test verts == cart.([(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)])

    # random points in interior do not affect result
    p1 = cart.([(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)])
    p2 = cart.([0.5 .* (rand(), rand()) .+ 0.5 for _ in 1:10])
    pts = [p1; p2]
    chul = hull(pts, method)
    verts = vertices(chul)
    @test verts == cart.([(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)])

    pts = cart.([
      (0, 5),
      (1, 5),
      (1, 4),
      (2, 4),
      (2, 3),
      (3, 3),
      (4, 3),
      (5, 3),
      (5, 4),
      (6, 4),
      (6, 5),
      (7, 5),
      (7, 6),
      (7, 7),
      (6, 7),
      (6, 8),
      (5, 8),
      (5, 9),
      (4, 9),
      (3, 9),
      (2, 9),
      (2, 8),
      (1, 8),
      (1, 7),
      (0, 7),
      (0, 6)
    ])
    chul = hull(pts, method)
    @test nvertices(chul) < length(pts)

    poly = readpoly(T, joinpath(datadir, "hull.line"))
    pts = vertices(poly)
    chul = hull(pts, method)
    @test nvertices(chul) < length(pts)

    if method == GrahamScan()
      # simplifying rectangular hull / triangular
      points = [cart(i - 1, j - 1) for i in 1:11 for j in 1:11]
      chull = hull(points, method)
      @test vertices(chull) == [cart(0, 0), cart(10, 0), cart(10, 10), cart(0, 10)]
      for _ in 1:100 # test presence of interior points doesn't affect the result
        push!(points, cart(10 * rand(), 10 * rand()))
      end
      chull = hull(points, method)
      @test vertices(chull) == [cart(0, 0), cart(10, 0), cart(10, 10), cart(0, 10)]

      points = [cart(-1, 0), cart(0, 0), cart(1, 0), cart(0, 2)]
      chull = hull(points, method)
      @test vertices(chull) == [cart(-1, 0), cart(1, 0), cart(0, 2)]

      # degenerate cases
      points = [cart(0, 0), cart(1, 0), cart(2, 0)]
      chull = hull(points, method)
      @test vertices(chull) == [cart(0, 0), cart(2, 0)]

      points = [cart(0, 0), cart(1, 0), cart(2, 0), cart(10, 0), cart(100, 0)]
      chull = hull(points, method)
      @test vertices(chull) == [cart(0, 0), cart(100, 0)]

      # partially collinear
      points = [
        cart(2, 0),
        cart(4, 0),
        cart(6, 0),
        cart(10, 0),
        cart(12, 1),
        cart(14, 3),
        cart(14, 6),
        cart(14, 9),
        cart(13, 10),
        cart(11, 11),
        cart(8, 12),
        cart(3, 11),
        cart(0, 8),
        cart(0, 7),
        cart(0, 6),
        cart(0, 5),
        cart(0, 4),
        cart(0, 3),
        cart(0, 2),
        cart(1, 0)
      ]
      chull = hull(points, method)
      truth = [
        cart(0, 2),
        cart(1, 0),
        cart(10, 0),
        cart(12, 1),
        cart(14, 3),
        cart(14, 9),
        cart(13, 10),
        cart(11, 11),
        cart(8, 12),
        cart(3, 11),
        cart(0, 8),
        cart(0, 3)
      ]
      @test vertices(chull) == truth
      push!(points, cart(4, 8), cart(2, 6), cart(6, 2), cart(10, 8), cart(8, 8), cart(10, 6))
      chull = hull(points, method)
      @test vertices(chull) == truth

      # https://github.com/JuliaGeometry/Meshes.jl/issues/1211
      data = readdlm(joinpath(datadir, "issue1211.dat"))
      points = cart.(data[:, 1], data[:, 2])
      chull = hull(points, method)
      @test area(chull) ≈ T(0.0015160200648848573)u"m^2"
    end
  end
end

@testitem "Convex hulls" setup = [Setup] begin
  @test convexhull(cart(0, 0)) == cart(0, 0)

  @test convexhull(Box(cart(0, 0), cart(1, 1))) == Box(cart(0, 0), cart(1, 1))

  @test convexhull(Ball(cart(0, 0), T(1))) == Ball(cart(0, 0), T(1))
  @test convexhull(Ball(cart(1, 1), T(1))) == Ball(cart(1, 1), T(1))

  @test convexhull(Sphere(cart(0, 0), T(1))) == Ball(cart(0, 0), T(1))
  @test convexhull(Sphere(cart(1, 1), T(1))) == Ball(cart(1, 1), T(1))

  b1 = Box(cart(0, 0), cart(1, 1))
  b2 = Box(cart(-1, -1), cart(0.5, 0.5))
  @test convexhull(Multi([b1, b2])) == PolyArea(cart.([(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)]))
  @test convexhull(GeometrySet([b1, b2])) == PolyArea(cart.([(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)]))

  b1 = Ball(cart(0, 0), T(1))
  b2 = Box(cart(-1, -1), cart(0, 0))
  h = convexhull(Multi([b1, b2]))
  @test cart(-0.8, -0.8) ∈ h
  @test cart(0.2, 0.2) ∈ h
end

@testitem "MoreiraSantosMarch" setup = [Setup] begin
  # constructor validation: k must be an integer > 2
  @test_throws AssertionError MoreiraSantosMarch(2)

  # basic random points
  pts = [cart(rand(T), rand(T)) for _ in 1:10]
  chul = hull(pts, MoreiraSantosMarch(3))
  @test all(pts .∈ Ref(chul))

  # corner cases bypass k entirely
  pt = cart.([(0, 0)])
  @test hull(pt, MoreiraSantosMarch(3)) == cart(0, 0)
  line = cart.([(0, 0), (1, 0)])
  @test hull(line, MoreiraSantosMarch(3)) == Segment(cart(0, 0), cart(1, 0))
  triangle = cart.([(1, 0), (0, 0), (0, 1)])
  chul = hull(triangle, MoreiraSantosMarch(3))
  @test Set(vertices(chul)) == Set(triangle)

  # hull-level k < n validation still applies once n > 3
  pts = cart.([(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)])
  @test_throws AssertionError hull(pts, MoreiraSantosMarch(5))
  @test_throws AssertionError hull(pts, MoreiraSantosMarch(6))
  chul = hull(pts, MoreiraSantosMarch(4))
  @test issimple(chul) && all(pts .∈ Ref(chul))

  # fuzz test without a try/catch escape hatch: must never throw
  rng = StableRNG(123)
  for _ in 1:100, k in (3, 4, 5)
    rpts = [cart(rand(rng, T), rand(rng, T)) for _ in 1:10]
    local chul = hull(rpts, MoreiraSantosMarch(k))
    @test nvertices(chul) ≥ 3
    @test all(rpts .∈ Ref(chul))
  end

  # U-shaped point set with a notch between x=1 and x=3 above y=1
  pts = cart.([
    (0, 0),
    (1, 0),
    (2, 0),
    (3, 0),
    (4, 0),
    (4, 1),
    (4, 2),
    (4, 3),
    (4, 4),
    (3, 4),
    (3, 3),
    (3, 2),
    (3, 1),
    (2, 1),
    (1, 1),
    (1, 2),
    (1, 3),
    (1, 4),
    (0, 4),
    (0, 3),
    (0, 2),
    (0, 1)
  ])
  chul = hull(pts, MoreiraSantosMarch(3))
  @test issimple(chul) && nvertices(chul) ≥ 3
  @test all(pts .∈ Ref(chul))
  @test Set(vertices(chul)) == Set(pts)
  @test area(chul) ≈ T(10) * u"m^2"
  @test area(chul) < area(hull(pts, JarvisMarch()))

  # moreira self-intersection regression
  pts = cart.([
    (3.7, 12.9),
    (5.9, 12.9),
    (9.3, 12.9),
    (10.4, 11.8),
    (1.5, 10.7),
    (7.9, 10.4),
    (0.4, 8.4),
    (3.0, 8.2),
    (5.7, 8.2),
    (4.4, 6.0),
    (0.4, 5.1),
    (1.5, 2.9),
    (7.0, 0.6),
    (5.7, 3.8),
    (9.3, 2.9),
    (4.8, 1.7)
  ])
  chul = hull(pts, MoreiraSantosMarch(3))
  @test issimple(chul) && nvertices(chul) ≥ 3
  @test all(pts .∈ Ref(chul))

  # true concavity regression against real data
  pts = vertices(readpoly(T, joinpath(datadir, "hull.line")))
  chul3 = hull(pts, MoreiraSantosMarch(3))
  chuln = hull(pts, MoreiraSantosMarch(length(pts) - 1))
  @test issimple(chul3) && all(pts .∈ Ref(chul3))
  @test issimple(chuln) && all(pts .∈ Ref(chuln))
  @test nvertices(chuln) ≤ nvertices(chul3)
end
