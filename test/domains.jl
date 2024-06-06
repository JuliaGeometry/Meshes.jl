@testset "Domain" begin
  # basic properties
  dom = DummyDomain(point(0, 0))
  @test embeddim(dom) == 2
  @test Meshes.crs(dom) <: Cartesian{NoDatum}
  @test Meshes.lentype(dom) == ℳ
  @test !isparametrized(dom)

  # indexable/iterable interface
  dom = DummyDomain(point(0, 0))
  @test dom[begin] == Ball(point(1, 1), T(1))
  @test dom[end] == Ball(point(3, 3), T(1))
  @test eltype(dom) <: Ball{2}
  @test length(dom) == 3
  @test keys(dom) == 1:3
  @test collect(dom) == [Ball(point(i, i), T(1)) for i in 1:3]
  @test dom[1:2] == [Ball(point(i, i), T(1)) for i in 1:2]

  # coordinates of centroids
  dom = DummyDomain(point(1, 1))
  pts = centroid.(Ref(dom), 1:3)
  @test pts == point.([(2, 2), (3, 3), (4, 4)])

  # concatenation
  dom1 = DummyDomain(point(0, 0))
  dom2 = DummyDomain(point(3, 3))
  dom3 = PointSet(randpoint2(3))
  @test vcat(dom1, dom2) == GeometrySet([collect(dom1); collect(dom2)])
  @test vcat(dom2, dom3) == GeometrySet([collect(dom2); collect(dom3)])
  @test vcat(dom3, dom1) == GeometrySet([collect(dom3); collect(dom1)])
  @test vcat(dom1, dom2, dom3) == GeometrySet([collect(dom1); collect(dom2); collect(dom3)])

  # datum propagation
  c = Cartesian{WGS84Latest}(T(1), T(1))
  dom = DummyDomain(Point(c))
  @test datum(Meshes.crs(centroid(dom))) === WGS84Latest

  dom = DummyDomain(point(0, 0))
  @test sprint(show, dom) == "3 DummyDomain"
  @test sprint(show, MIME"text/plain"(), dom) == """
  3 DummyDomain
  ├─ Ball(center: (x: 1.0 m, y: 1.0 m), radius: 1.0 m)
  ├─ Ball(center: (x: 2.0 m, y: 2.0 m), radius: 1.0 m)
  └─ Ball(center: (x: 3.0 m, y: 3.0 m), radius: 1.0 m)"""
end
