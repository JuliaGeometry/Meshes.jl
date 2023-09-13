# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------
# DOMAIN VIEWS
# -------------

"""
    SubDomain(domain, indices)

A partial view of a `domain` containing only the elements at `indices`.
"""
struct SubDomain{Dim,T,D<:Domain{Dim,T},I} <: Domain{Dim,T}
  domain::D
  inds::I
end

# specialize constructor to avoid infinite loops
SubDomain(v::SubDomain, inds) = SubDomain(getfield(v, :domain), getfield(v, :inds)[inds])

# -----------------
# DOMAIN INTERFACE
# -----------------

element(v::SubDomain, ind::Int) = element(v.domain, v.inds[ind])

nelements(v::SubDomain) = length(v.inds)

centroid(v::SubDomain, ind::Int) = centroid(v.domain, v.inds[ind])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::SubDomain)
  domain = getfield(v, :domain)
  nelms = length(getfield(v, :inds))
  print(io, "$nelms View{$domain}")
end
