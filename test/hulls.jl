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
    @test chul == Segment(P2(1, 0), P2(0, 1))
    pset = P2[(1, 0), (0, 0), (0, 1)]
    chul = hull(pset, GrahamScan())
    @test vertices(boundary(chul)) == P2[(0, 0), (1, 0), (0, 1)]

    # original point set is already in hull
    pset = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)]
    chul = hull(pset, GrahamScan())
    verts = vertices(boundary(chul))
    @test verts == P2[(0.5, -1), (1, 0), (1, 1), (0, 1), (0, 0)]

    # random points in interior do not affect result
    p1 = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0.5, -1)]
    p2 = P2[0.5 .* (rand(), rand()) .+ 0.5 for _ in 1:10]
    pset = [p1; p2]
    chul = hull(pset, GrahamScan())
    verts = vertices(boundary(chul))
    @test verts == P2[(0.5, -1), (1, 0), (1, 1), (0, 1), (0, 0)]
  end

  @testset "Miscelanneous" begin
    @test hull(P2(0, 0)) == Box(P2(0, 0), P2(0, 0))

    @test hull(Box(P2(0, 0), P2(1, 1))) == Box(P2(0, 0), P2(1, 1))

    @test hull(Ball(P2(0, 0), T(1))) == Ball(P2(0, 0), T(1))
    @test hull(Ball(P2(1, 1), T(1))) == Ball(P2(1, 1), T(1))

    @test hull(Sphere(P2(0, 0), T(1))) == Ball(P2(0, 0), T(1))
    @test hull(Sphere(P2(1, 1), T(1))) == Ball(P2(1, 1), T(1))

    b1 = Box(P2(0, 0), P2(1, 1))
    b2 = Box(P2(-1, -1), P2(0.5, 0.5))
    @test hull(Multi([b1, b2])) == PolyArea(P2[(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)])
    @test hull(Collection([b1, b2])) == PolyArea(P2[(-1, -1), (0.5, -1), (1, 0), (1, 1), (0, 1), (-1, 0.5)])
  end
end
