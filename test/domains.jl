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

  dom = DummyDomain(P2(0, 0))
  @test sprint(show, dom) == "3 DummyDomain{2,$T}"
  if T == Float32
    @test sprint(show, MIME"text/plain"(), dom) == """
    3 DummyDomain{2,Float32}
    ├─ Ball(center = (1.0f0, 1.0f0), radius = 1.0)
    ├─ Ball(center = (2.0f0, 2.0f0), radius = 1.0)
    └─ Ball(center = (3.0f0, 3.0f0), radius = 1.0)"""
  elseif T == Float64
    @test sprint(show, MIME"text/plain"(), dom) == """
    3 DummyDomain{2,Float64}
    ├─ Ball(center = (1.0, 1.0), radius = 1.0)
    ├─ Ball(center = (2.0, 2.0), radius = 1.0)
    └─ Ball(center = (3.0, 3.0), radius = 1.0)"""
  end
end
