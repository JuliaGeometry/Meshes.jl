@testitem "Hulls" setup = [Setup] begin
  for method in [GrahamScan(), JarvisMarch()]
    # basic test
    pts = randpoint2(100)
    chul = hull(pts, method)
    @test all(pts .∈ Ref(chul))

    # duplicated points
    pts = [randpoint2(100); randpoint2(100)]
    chul = hull(pts, method)
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

    pts =
      cart.([
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
      @test vertices(chull) == (cart(0, 0), cart(2, 0))

      points = [cart(0, 0), cart(1, 0), cart(2, 0), cart(10, 0), cart(100, 0)]
      chull = hull(points, method)
      @test vertices(chull) == (cart(0, 0), cart(100, 0))

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
        cart(0, 8)
      ]
      @test vertices(chull) == truth
      push!(points, cart(4, 8), cart(2, 6), cart(6, 2), cart(10, 8), cart(8, 8), cart(10, 6))
      chull = hull(points, method)
      @test vertices(chull) == truth
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
