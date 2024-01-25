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
