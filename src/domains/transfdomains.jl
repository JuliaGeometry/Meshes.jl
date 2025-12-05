# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TransformedDomain(domain, transform)

Lazy representation of a geometric `transform` applied to a `domain`.
"""
struct TransformedDomain{M<:Manifold,C<:CRS,D<:Domain,T<:Transform} <: Domain{M,C}
  domain::D
  transform::T

  function TransformedDomain{M,C}(domain::D, transform::T) where {M<:Manifold,C<:CRS,D<:Domain,T<:Transform}
    new{M,C,D,T}(domain, transform)
  end
end

function TransformedDomain(d::Domain, t::Transform)
  g = t(first(d))
  TransformedDomain{manifold(g),crs(g)}(d, t)
end

# specialize constructor to avoid deep structures
TransformedDomain(d::TransformedDomain, t::Transform) = TransformedDomain(d.domain, d.transform → t)

Base.parent(d::TransformedDomain) = d.domain

transform(d::TransformedDomain) = d.transform

# domain interface
element(d::TransformedDomain, ind::Int) = d.transform(element(d.domain, ind))

nelements(d::TransformedDomain) = nelements(d.domain)
