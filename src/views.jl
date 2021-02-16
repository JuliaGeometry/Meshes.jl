# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DomainView(domain, indices)

A partial view of a `domain` containing only the elements at `indices`.
"""
struct DomainView{Dim,T,D<:Domain{Dim,T},I} <: Domain{Dim,T}
  domain::D
  inds::I
end

# convenience functions
Base.view(domain::Domain, inds) = DomainView(domain, inds)

# -----------------
# DOMAIN INTERFACE
# -----------------

Base.getindex(v::DomainView, ind::Int) =
  getindex(v.domain, v.inds[ind])

nelements(v::DomainView) = length(v.inds)

coordinates!(buff, v::DomainView, ind::Int) = 
  coordinates!(buff, v.domain, v.inds[ind])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::DomainView)
  domain  = v.domain
  nelms = length(v.inds)
  print(io, "$nelms View{$domain}")
end
