@testset "Hulls" begin
  @testset "GrahamScan" begin
    # basic test
    pset = rand(P2, 100)
    chul = hull(pset, GrahamScan())
    @test all(pset .∈ Ref(chul))

    # duplicated points
    pset = [rand(P2, 100); rand(P2, 100)]
    chul = hull(pset, GrahamScan())
    @test all(pset .∈ Ref(chul))

    # corner cases
    pset = P2[(0, 0)]
    chul = hull(pset, GrahamScan())
    @test chul == P2(0, 0)
    pset = P2[(0, 1), (1, 0)]
    chul = hull(pset, GrahamScan())
    @test chul == Segment(P2(0, 1), P2(1, 0))
    pset = P2[(1, 0), (0, 0), (0, 1)]
    chul = hull(pset, GrahamScan())
    @test vertices(chul) == P2[(0, 0), (1, 0), (0, 1)]

    # original point set is already in hull
    pset = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)]
    chul = hull(pset, GrahamScan())
    verts = vertices(chul)
    @test verts == P2[(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)]

    # random points in interior do not affect result
    p1 = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)]
    p2 = P2[0.5 .* (rand(), rand()) .+ 0.5 for _ in 1:10]
    pset = [p1; p2]
    chul = hull(pset, GrahamScan())
    verts = vertices(chul)
    @test verts == P2[(0, 0), (0.5, -1), (1, 0), (1, 1), (0, 1)]
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
