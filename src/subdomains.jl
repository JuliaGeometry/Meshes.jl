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
struct SubDomain{Dim,T,D<:Domain{Dim,T},I<:AbstractVector{Int}} <: Domain{Dim,T}
  domain::D
  inds::I
end

# specialize constructor to avoid infinite loops
SubDomain(v::SubDomain, inds::AbstractVector{Int}) = SubDomain(v.domain, v.inds[inds])

# -----------------
# DOMAIN INTERFACE
# -----------------

element(v::SubDomain, ind::Int) = element(v.domain, v.inds[ind])

nelements(v::SubDomain) = length(v.inds)

centroid(v::SubDomain, ind::Int) = centroid(v.domain, v.inds[ind])

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, v::SubDomain{Dim,T}) where {Dim,T}
  domain = getfield(v, :domain)
  name = prettyname(domain)
  inds = getfield(v, :inds)
  nelms = length(inds)
  print(io, "$nelms view(::$name{$Dim,$T}, ")
  printinds(io, inds)
  print(io, ")")
end
