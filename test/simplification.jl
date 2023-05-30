@testset "Simplification" begin
  @testset "DouglasPeucker" begin
    c = Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])
    s1 = simplify(c, DouglasPeucker(T(0.1)))
    s2 = simplify(c, DouglasPeucker(T(0.5)))
    @test s1 == Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])
    @test s2 == Ring(P2[(0, 0), (1.5, 0.5), (0, 1)])

    p = PolyArea(Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)]))
    s1 = simplify(p, DouglasPeucker(T(0.5)))
    @test s1 == PolyArea(Ring(P2[(0, 0), (1.5, 0.5), (0, 1)]))
    m = Multi([p, p])
    s2 = simplify(m, DouglasPeucker(T(0.5)))
    @test s2 == Multi([s1, s1])
    d = GeometrySet([p, p])
    s3 = simplify(d, DouglasPeucker(T(0.5)))
    @test s3 == GeometrySet([s1, s1])

    # perform binary search for ϵ tolerance
    c = Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])
    s1 = simplify(c, DouglasPeucker(T(0.1)))
    s2 = simplify(c, DouglasPeucker(max=6))
    @test s1 == s2
    s1 = simplify(c, DouglasPeucker(T(0.5)))
    s2 = simplify(c, DouglasPeucker(max=4))
    @test s1 == s2
  end

  @testset "Selinger" begin
    c = Ring(P2[(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (0, 2), (0, 1)])
    s1 = simplify(c, Selinger(0.1))
    s2 = simplify(c, Selinger(0.5))
    @test s1 == Ring(P2[(1, 0), (1, 1), (2, 1), (2, 2), (0, 2), (0, 0)])
    @test s2 == Ring(P2[(1, 0), (2, 2), (0, 2), (0, 0)])
  end

  @testset "Utilities" begin
    # decimate is a helper function to simplify
    # geometries with an appropriate method
    b = Box(P2(0, 0), P2(1, 1))
    s = decimate(b, 1.0)
    @test s isa Polygon
    @test nvertices(s) == 3
    @test boundary(s) == Ring(P2[(0, 0), (1, 0), (0, 1)])

    c = Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])
    s1 = decimate(c, T(0.1))
    s2 = decimate(c, T(0.5))
    @test s1 == Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])
    @test s2 == Ring(P2[(0, 0), (1.5, 0.5), (0, 1)])

    c = Ring(P2[(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])
    s1 = decimate(c, T(0.1))
    s2 = decimate(c, max=6)
    @test s1 == s2
    s1 = decimate(c, T(0.5))
    s2 = decimate(c, max=4)
    @test s1 == s2
  end
end
