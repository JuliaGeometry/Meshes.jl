@testitem "Domain" setup = [Setup] begin
  # basic properties
  dom = DummyDomain(cart(0, 0))
  @test embeddim(dom) == 2
  @test crs(dom) <: Cartesian{NoDatum}
  @test Meshes.lentype(dom) == ℳ
  @test !isparametrized(dom)

  # indexable/iterable interface
  dom = DummyDomain(cart(0, 0))
  @test dom[begin] == Ball(cart(1, 1), T(1))
  @test dom[end] == Ball(cart(3, 3), T(1))
  @test eltype(dom) <: Ball
  @test length(dom) == 3
  @test keys(dom) == 1:3
  @test collect(dom) == [Ball(cart(i, i), T(1)) for i in 1:3]
  @test dom[1:2] == [Ball(cart(i, i), T(1)) for i in 1:2]

  # coordinates of centroids
  dom = DummyDomain(cart(1, 1))
  pts = centroid.(Ref(dom), 1:3)
  @test pts == cart.([(2, 2), (3, 3), (4, 4)])

  # concatenation
  dom1 = DummyDomain(cart(0, 0))
  dom2 = DummyDomain(cart(3, 3))
  dom3 = PointSet([cart(1, 1), cart(2, 2), cart(3, 3)])
  @test vcat(dom1, dom2) == GeometrySet([collect(dom1); collect(dom2)])
  @test vcat(dom2, dom3) == GeometrySet([collect(dom2); collect(dom3)])
  @test vcat(dom3, dom1) == GeometrySet([collect(dom3); collect(dom1)])
  @test vcat(dom1, dom2, dom3) == GeometrySet([collect(dom1); collect(dom2); collect(dom3)])

  # CRS propagation
  dom = DummyDomain(merc(1, 1))
  @test crs(centroid(dom)) === crs(dom)

  dom = DummyDomain(cart(0, 0))
  @test sprint(show, dom) == "3 DummyDomain"
  @test sprint(show, MIME"text/plain"(), dom) == """
  3 DummyDomain
  ├─ Ball(center: (x: 1.0 m, y: 1.0 m), radius: 1.0 m)
  ├─ Ball(center: (x: 2.0 m, y: 2.0 m), radius: 1.0 m)
  └─ Ball(center: (x: 3.0 m, y: 3.0 m), radius: 1.0 m)"""
end
