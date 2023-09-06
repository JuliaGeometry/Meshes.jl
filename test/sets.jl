@testset "Sets" begin
  @testset "GeometrySet" begin
    s = Segment(P2(0, 0), P2(1, 1))
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    p = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    gset = GeometrySet([s, t, p])
    @test [centroid(gset, i) for i in 1:3] == P2[(1 / 2, 1 / 2), (1 / 3, 1 / 3), (1 / 2, 1 / 2)]

    s = Segment(P2(0, 0), P2(1, 1))
    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    geoms = [s, t]
    gset1 = GeometrySet(geoms)
    gset2 = GeometrySet(g for g in geoms)
    @test gset1 == gset2
    @test parent(gset1) === geoms
  end

  @testset "PointSet" begin
    pset = PointSet(rand(P1, 100))
    @test embeddim(pset) == 1
    @test coordtype(pset) == T
    @test nelements(pset) == 100
    @test eltype(pset) <: P1

    pset = PointSet(rand(P2, 100))
    @test embeddim(pset) == 2
    @test coordtype(pset) == T
    @test nelements(pset) == 100
    @test eltype(pset) <: P2

    pset = PointSet(rand(P3, 100))
    @test embeddim(pset) == 3
    @test coordtype(pset) == T
    @test nelements(pset) == 100
    @test eltype(pset) <: P3

    pset1 = PointSet([P3(1, 2, 3), P3(4, 5, 6)])
    pset2 = PointSet(P3(1, 2, 3), P3(4, 5, 6))
    pset3 = PointSet([T.((1, 2, 3)), T.((4, 5, 6))])
    pset4 = PointSet(T.((1, 2, 3)), T.((4, 5, 6)))
    pset5 = PointSet(T[1 4; 2 5; 3 6])
    @test pset1 == pset2 == pset3 == pset4 == pset5
    for pset in [pset1, pset2, pset3, pset4, pset5]
      @test embeddim(pset) == 3
      @test coordtype(pset) == T
      @test nelements(pset) == 2
      @test pset[1] == P3(1, 2, 3)
      @test pset[2] == P3(4, 5, 6)
    end

    pset = PointSet(P2[(0, 0), (1, 0), (0, 1)])
    @test centroid(pset) == P2(1 / 3, 1 / 3)

    pset = PointSet(P2[(1, 0), (0, 1)])
    @test nelements(pset) == 2
    @test centroid(pset, 1) == P2(1, 0)
    @test centroid(pset, 2) == P2(0, 1)

    pset = PointSet(P2[(0, 0), (1, 0), (0, 1)])
    @test measure(pset) == zero(T)

    # constructor with iterator
    points = P2[(1, 0), (0, 1)]
    pset1 = PointSet(points)
    pset2 = PointSet(p for p in points)
    @test pset1 == pset2

    pset = PointSet(P2[(1, 0), (0, 1)])
    @test sprint(show, pset) == "2 PointSet{2,$T}"
    if T == Float32
      @test sprint(show, MIME"text/plain"(), pset) == """
      2 PointSet{2,Float32}
      ├─ Point(1.0f0, 0.0f0)
      └─ Point(0.0f0, 1.0f0)"""
    elseif T == Float64
      @test sprint(show, MIME"text/plain"(), pset) == """
      2 PointSet{2,Float64}
      ├─ Point(1.0, 0.0)
      └─ Point(0.0, 1.0)"""
    end
  end
end
