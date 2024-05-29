@testset "Hulls" begin
  @testset "Basic" begin
    for method in [GrahamScan(), JarvisMarch()]
      # basic test
      pts = rand(Point{2}, 100)
      chul = hull(pts, method)
      @test all(pts .∈ Ref(chul))

      # duplicated points
      pts = [rand(Point{2}, 100); rand(Point{2}, 100)]
      chul = hull(pts, method)
      @test all(pts .∈ Ref(chul))

      # corner cases
      pts = point.([(0, 0)])
      chul = hull(pts, method)
      @test chul == point(0, 0)
      pts = point.([(0, 1), (1, 0)])
      chul = hull(pts, method)
      @test chul == Segment(point(0, 1), point(1, 0))
      pts = point.([(1, 0), (0, 0), (0, 1)])
      chul = hull(pts, method)
      @test vertices(chul) == point.([(0, 0), (1, 0), (0, 1)])

      # original point set is already in hull
      pts = point.([(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)])
      chul = hull(pts, method)
      verts = vertices(chul)
      @test verts == point.([(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)])

      # random points in interior do not affect result
      p1 = point.([(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)])
      p2 = point.([0.5 .* (rand(), rand()) .+ 0.5 for _ in 1:10])
      pts = [p1; p2]
      chul = hull(pts, method)
      verts = vertices(chul)
      @test verts == point.([(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)])

      pts =
        point.([
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
        points = [point(i - 1, j - 1) for i in 1:11 for j in 1:11]
        chull = hull(points, method)
        @test vertices(chull) == [point(0, 0), point(10, 0), point(10, 10), point(0, 10)]
        for _ in 1:100 # test presence of interior points doesn't affect the result 
          push!(points, point(10 * rand(), 10 * rand()))
        end
        chull = hull(points, method)
        @test vertices(chull) == [point(0, 0), point(10, 0), point(10, 10), point(0, 10)]

        points = [point(-1, 0), point(0, 0), point(1, 0), point(0, 2)]
        chull = hull(points, method)
        @test vertices(chull) == [point(-1, 0), point(1, 0), point(0, 2)]

        # degenerate cases
        points = [point(0, 0), point(1, 0), point(2, 0)]
        chull = hull(points, method)
        @test vertices(chull) == (point(0, 0), point(2, 0))

        points = [point(0, 0), point(1, 0), point(2, 0), point(10, 0), point(100, 0)]
        chull = hull(points, method)
        @test vertices(chull) == (point(0, 0), point(100, 0))

        # partially collinear 
        points = [
          point(2, 0),
          point(4, 0),
          point(6, 0),
          point(10, 0),
          point(12, 1),
          point(14, 3),
          point(14, 6),
          point(14, 9),
          point(13, 10),
          point(11, 11),
          point(8, 12),
          point(3, 11),
          point(0, 8),
          point(0, 7),
          point(0, 6),
          point(0, 5),
          point(0, 4),
          point(0, 3),
          point(0, 2),
          point(1, 0)
        ]
        chull = hull(points, method)
        truth = [
          point(0, 2),
          point(1, 0),
          point(10, 0),
          point(12, 1),
          point(14, 3),
          point(14, 9),
          point(13, 10),
          point(11, 11),
          point(8, 12),
          point(3, 11),
          point(0, 8)
        ]
        @test vertices(chull) == truth
        push!(points, point(4, 8), point(2, 6), point(6, 2), point(10, 8), point(8, 8), point(10, 6))
        chull = hull(points, method)
        @test vertices(chull) == truth
      end
    end
  end

  @testset "convexhull" begin
    @test convexhull(point(0, 0)) == point(0, 0)

    @test convexhull(Box(point(0, 0), point(1, 1))) == Box(point(0, 0), point(1, 1))

    @test convexhull(Ball(point(0, 0), T(1))) == Ball(point(0, 0), T(1))
    @test convexhull(Ball(point(1, 1), T(1))) == Ball(point(1, 1), T(1))

    @test convexhull(Sphere(point(0, 0), T(1))) == Ball(point(0, 0), T(1))
    @test convexhull(Sphere(point(1, 1), T(1))) == Ball(point(1, 1), T(1))

    b1 = Box(point(0, 0), point(1, 1))
    b2 = Box(point(-1, -1), point(0.5, 0.5))
    @test convexhull(Multi([b1, b2])) == PolyArea(point.([(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)]))
    @test convexhull(GeometrySet([b1, b2])) ==
          PolyArea(point.([(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)]))

    b1 = Ball(point(0, 0), T(1))
    b2 = Box(point(-1, -1), point(0, 0))
    h = convexhull(Multi([b1, b2]))
    @test point(-0.8, -0.8) ∈ h
    @test point(0.2, 0.2) ∈ h
  end
end
