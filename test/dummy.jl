# dummy type implementing the Domain trait
struct DummyDomain{Dim,P<:Point{Dim}} <: Domain{Dim}
  origin::P
end

Meshes.lentype(::Type{<:DummyDomain{Dim,P}}) where {Dim,P} = Meshes.lentype(P)

function Meshes.element(domain::DummyDomain{Dim}, ind::Int) where {Dim}
  ℒ = Meshes.lentype(domain)
  T = Unitful.numtype(ℒ)
  c = domain.origin + Vec(ntuple(i -> T(ind) * unit(ℒ), Dim))
  r = oneunit(ℒ)
  Ball(c, r)
end

Meshes.nelements(d::DummyDomain) = 3
