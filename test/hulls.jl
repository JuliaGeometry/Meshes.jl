@testset "Hulls" begin
  @testset "GrahamScan" begin
    # basic test
    pset = PointSet(rand(P2, 100))
    chul = hull(pset, GrahamScan())
    @test all(pset .∈ Ref(chul))

    # duplicated points
    pset = PointSet([rand(P2, 100); rand(P2, 100)])
    chul = hull(pset, GrahamScan())
    @test all(pset .∈ Ref(chul))

    # corner cases
    pset = PointSet(P2[(0,0)])
    chul = hull(pset, GrahamScan())
    @test chul == P2(0,0)
    pset = PointSet(P2[(0,1),(1,0)])
    chul = hull(pset, GrahamScan())
    @test chul == Segment(P2(1,0), P2(0,1))
    pset = PointSet(P2[(1,0),(0,0),(0,1)])
    chul = hull(pset, GrahamScan())
    @test vertices(boundary(chul)) == P2[(0,0),(1,0),(0,1)]

    # original point set is already in hull
    pset  = PointSet(P2[(0,0),(1,0),(1,1),(0,1),(0.5,-1)])
    chul  = hull(pset, GrahamScan())
    verts = vertices(boundary(chul))
    @test verts == P2[(0.5,-1),(1,0),(1,1),(0,1),(0,0)]

    # random points in interior do not affect result
    p1 = P2[(0,0),(1,0),(1,1),(0,1),(0.5,-1)]
    p2 = P2[0.5.*(rand(), rand()) .+ 0.5 for _ in 1:10]
    pset  = PointSet([p1; p2])
    chul  = hull(pset, GrahamScan())
    verts = vertices(boundary(chul))
    @test verts == P2[(0.5,-1),(1,0),(1,1),(0,1),(0,0)]
  end
end