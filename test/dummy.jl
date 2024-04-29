# dummy type implementing the Domain trait
struct DummyDomain{Dim,P<:Point{Dim}} <: Domain{Dim}
  origin::P
end
function Meshes.element(domain::DummyDomain{Dim}, ind::Int) where {Dim}
  c = domain.origin + Vec(ntuple(i -> ind, Dim))
  r = one(T)
  Ball(c, r)
end
Meshes.nelements(d::DummyDomain) = 3
