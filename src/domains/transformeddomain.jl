# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TranformedDomain(domain, transform)

TODO
"""
struct TranformedDomain{M<:Manifold,C<:CRS,D<:Domain,T<:Transform} <: Domain{M,C}
  domain::D
  transform::T

  function TranformedDomain{M,C}(domain::D, transform::T) where {M<:Manifold,C<:CRS,D<:Domain,T<:Transform}
    new{M,C,D,T}(domain, transform)
  end
end

function TranformedDomain(d::Domain, t::Transform)
  g = t(first(d))
  TransformedMesh{manifold(g),crs(g)}(d, t)
end

# specialize constructor to avoid deep structures
TranformedDomain(d::TranformedDomain, t::Transform) = TransformedMesh(d.domain, d.transform → t)

Base.parent(d::TranformedDomain) = d.domain

transform(d::TranformedDomain) = d.transform

# domain interface
element(d::TranformedDomain, ind::Int) = d.transform(element(d.domain, ind))

nelements(d::TranformedDomain) = nelements(d.domain)
