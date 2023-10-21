@testset "Hulls" begin
  @testset "Basic" begin
    for method in [GrahamScan(), JarvisMarch()]
      # basic test
      pts = rand(P2, 100)
      chul = hull(pts, method)
      @test all(pts .∈ Ref(chul))

      # duplicated points
      pts = [rand(P2, 100); rand(P2, 100)]
      chul = hull(pts, method)
      @test all(pts .∈ Ref(chul))

      # corner cases
      pts = P2[(0, 0)]
      chul = hull(pts, method)
      @test chul == P2(0, 0)
      pts = P2[(0, 1), (1, 0)]
      chul = hull(pts, method)
      @test chul == Segment(P2(0, 1), P2(1, 0))
      pts = P2[(1, 0), (0, 0), (0, 1)]
      chul = hull(pts, method)
      @test vertices(chul) == P2[(0, 0), (1, 0), (0, 1)]

      # original point set is already in hull
      pts = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)]
      chul = hull(pts, method)
      verts = vertices(chul)
      @test verts == P2[(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)]

      # random points in interior do not affect result
      p1 = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)]
      p2 = P2[0.5 .* (rand(), rand()) .+ 0.5 for _ in 1:10]
      pts = [p1; p2]
      chul = hull(pts, method)
      verts = vertices(chul)
      @test verts == P2[(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)]

      pts = P2[
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
      ]
      chul = hull(pts, method)
      @test nvertices(chul) < length(pts)

      poly = readpoly(T, joinpath(datadir, "hull.line"))
      pts = vertices(poly)
      chul = hull(pts, method)
      @test nvertices(chul) < length(pts)

      if method == GrahamScan() # simplification not yet implemented for JarvisMarch()
        # simplifying rectangular hull / triangular
        points = [P2(i - 1, j - 1) for i in 1:11, j in 1:11] |> vec
        chull = hull(points, method)
        @test vertices(chull) == [P2(0, 0), P2(10, 0), P2(10, 10), P2(0, 10)]
        for _ in 1:100 # test presence of interior points doesn't affect the result 
          push!(points, P2(10 * rand(), 10 * rand()))
        end
        chull = hull(points, method)
        @test vertices(chull) == [P2(0, 0), P2(10, 0), P2(10, 10), P2(0, 10)]

        points = [P2(-1, 0), P2(0, 0), P2(1, 0), P2(0, 2)]
        chull = hull(points, method)
        @test vertices(chull) == [P2(-1, 0), P2(1, 0), P2(0, 2)]

        # degenerate cases
        points = [P2(0, 0), P2(1, 0), P2(2, 0)]
        chull = hull(points, method)
        @test vertices(chull) == [P2(0, 0), P2(1, 0), P2(2, 0)] # the degeneracy fix of PolyArea inserts P2(1,0)

        points = [P2(0, 0), P2(1, 0), P2(2, 0), P2(10, 0), P2(100, 0)]
        chull = hull(points, method)
        @test vertices(chull) == [P2(0, 0), P2(50, 0), P2(100, 0)] # the degeneracy fix of PolyArea inserts P2(50,0)

        # partially collinear 
        points = [
          P2(2, 0),
          P2(4, 0),
          P2(6, 0),
          P2(10, 0),
          P2(12, 1),
          P2(14, 3),
          P2(14, 6),
          P2(14, 9),
          P2(13, 10),
          P2(11, 11),
          P2(8, 12),
          P2(3, 11),
          P2(0, 8),
          P2(0, 7),
          P2(0, 6),
          P2(0, 5),
          P2(0, 4),
          P2(0, 3),
          P2(0, 2),
          P2(1, 0)
        ]
        chull = hull(points, method)
        truth = [
          P2(0, 2),
          P2(1, 0),
          P2(10, 0),
          P2(12, 1),
          P2(14, 3),
          P2(14, 9),
          P2(13, 10),
          P2(11, 11),
          P2(8, 12),
          P2(3, 11),
          P2(0, 8)
        ]
        @test vertices(chull) == truth
        push!(points, P2(4, 8), P2(2, 6), P2(6, 2), P2(10, 8), P2(8, 8), P2(10, 6))
        chull = hull(points, method)
        @test vertices(chull) == truth
      end
    end
  end

  @testset "convexhull" begin
    @test convexhull(P2(0, 0)) == P2(0, 0)

    @test convexhull(Box(P2(0, 0), P2(1, 1))) == Box(P2(0, 0), P2(1, 1))

    @test convexhull(Ball(P2(0, 0), T(1))) == Ball(P2(0, 0), T(1))
    @test convexhull(Ball(P2(1, 1), T(1))) == Ball(P2(1, 1), T(1))

    @test convexhull(Sphere(P2(0, 0), T(1))) == Ball(P2(0, 0), T(1))
    @test convexhull(Sphere(P2(1, 1), T(1))) == Ball(P2(1, 1), T(1))

    b1 = Box(P2(0, 0), P2(1, 1))
    b2 = Box(P2(-1, -1), P2(0.5, 0.5))
    @test convexhull(Multi([b1, b2])) == PolyArea(P2[(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)])
    @test convexhull(GeometrySet([b1, b2])) == PolyArea(P2[(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)])

    b1 = Ball(P2(0, 0), T(1))
    b2 = Box(P2(-1, -1), P2(0, 0))
    h = convexhull(Multi([b1, b2]))
    @test P2(-0.8, -0.8) ∈ h
    @test P2(0.2, 0.2) ∈ h
  end
end