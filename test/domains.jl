@testset "Domain" begin
  # dummy type implementing the Domain trait
  struct DummyDomain{Dim,T} <: Domain{Dim,T}
    origin::Point{Dim,T}
  end
  function Meshes.element(domain::DummyDomain{Dim,T}, ind::Int) where {Dim,T}
    c = domain.origin + Vec(ntuple(i -> T(ind), Dim))
    r = one(T)
    Ball(c, r)
  end
  Meshes.nelements(d::DummyDomain) = 3

  # basic properties
  dom = DummyDomain(P2(0, 0))
  @test embeddim(dom) == 2
  @test coordtype(dom) == T
  @test !isparametrized(dom)

  # indexable/iterable interface
  dom = DummyDomain(P2(0, 0))
  @test dom[begin] == Ball(P2(1, 1), T(1))
  @test dom[end] == Ball(P2(3, 3), T(1))
  @test eltype(dom) <: Ball{2,T}
  @test length(dom) == 3
  @test keys(dom) == 1:3
  @test collect(dom) == [Ball(P2(i, i), T(1)) for i in 1:3]
  @test dom[1:2] == [Ball(P2(i, i), T(1)) for i in 1:2]

  # coordinates of centroids
  dom = DummyDomain(P2(1, 1))
  pts = centroid.(Ref(dom), 1:3)
  @test pts == P2[(2, 2), (3, 3), (4, 4)]

  # concatenation
  dom1 = DummyDomain(P2(0, 0))
  dom2 = DummyDomain(P2(3, 3))
  dom3 = PointSet(rand(P2, 3))
  @test vcat(dom1, dom2) == GeometrySet([collect(dom1); collect(dom2)])
  @test vcat(dom2, dom3) == GeometrySet([collect(dom2); collect(dom3)])
  @test vcat(dom3, dom1) == GeometrySet([collect(dom3); collect(dom1)])
  @test vcat(dom1, dom2, dom3) == GeometrySet([collect(dom1); collect(dom2); collect(dom3)])

  dom = DummyDomain(P2(0, 0))
  @test sprint(show, dom) == "3 DummyDomain{2,$T}"
  @test sprint(show, MIME"text/plain"(), dom) == """
  3 DummyDomain{2,$T}
  ├─ Ball(center: (1.0, 1.0), radius: 1.0)
  ├─ Ball(center: (2.0, 2.0), radius: 1.0)
  └─ Ball(center: (3.0, 3.0), radius: 1.0)"""
end
