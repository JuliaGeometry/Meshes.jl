# dummy type implementing the Domain trait
struct DummyDomain{Dim,C<:CRS} <: Domain{Dim,C}
  origin::Point{Dim,C}
end

function Meshes.element(domain::DummyDomain{Dim}, ind::Int) where {Dim}
  ℒ = Meshes.lentype(domain)
  T = Unitful.numtype(ℒ)
  c = domain.origin + Vec(ntuple(i -> T(ind) * unit(ℒ), Dim))
  r = oneunit(ℒ)
  Ball(c, r)
end

Meshes.nelements(d::DummyDomain) = 3
