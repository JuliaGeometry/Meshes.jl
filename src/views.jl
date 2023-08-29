# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# -------------
# DOMAIN VIEWS
# -------------

"""
    DomainView(domain, indices)

A partial view of a `domain` containing only the elements at `indices`.
"""
struct DomainView{Dim,T,D<:Domain{Dim,T},I} <: Domain{Dim,T}
  domain::D
  inds::I
end

# specialize constructor to avoid infinite loops
DomainView(v::DomainView, inds) = DomainView(getfield(v, :domain), getfield(v, :inds)[inds])

# -----------------
# DOMAIN INTERFACE
# -----------------

element(v::DomainView, ind::Int) = element(v.domain, v.inds[ind])

nelements(v::DomainView) = length(v.inds)

centroid(v::DomainView, ind::Int) = centroid(v.domain, v.inds[ind])

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, v::DomainView)
  domain = getfield(v, :domain)
  nelms = length(getfield(v, :inds))
  print(io, "$nelms View{$domain}")
end
