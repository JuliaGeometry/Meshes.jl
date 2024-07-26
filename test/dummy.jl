# dummy type implementing the Domain trait
struct DummyDomain{M<:Meshes.Manifold,C<:CRS} <: Domain{M,C}
  origin::Point{M,C}
end

function Meshes.element(domain::DummyDomain, ind::Int)
  ℒ = Meshes.lentype(domain)
  T = Unitful.numtype(ℒ)
  c = domain.origin + Vec(ntuple(i -> T(ind) * unit(ℒ), embeddim(domain)))
  r = oneunit(ℒ)
  Ball(c, r)
end

Meshes.nelements(d::DummyDomain) = 3
