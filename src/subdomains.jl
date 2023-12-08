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
SubDomain(d::SubDomain, inds::AbstractVector{Int}) = SubDomain(d.domain, d.inds[inds])

# -----------------
# DOMAIN INTERFACE
# -----------------

element(d::SubDomain, ind::Int) = element(d.domain, d.inds[ind])

nelements(d::SubDomain) = length(d.inds)

centroid(d::SubDomain, ind::Int) = centroid(d.domain, d.inds[ind])

# specialized for efficiency
function Base.vcat(d1::SubDomain, d2::SubDomain)
  if parent(d1) === parent(d2)
    dom = parent(d1)
    inds = vcat(parentindices(d1), parentindices(d2))
    SubDomain(dom, inds)
  else
    GeometrySet(vcat(collect(d1), collect(d2)))
  end
end

# -------------
# UNWRAP VIEWS
# -------------

"""
    parent(subdomain)

Returns the "parent domain" of a domain view.
"""
Base.parent(d::SubDomain) = d.domain

"""
    parentindices(subdomain)

Returns the indices used to create the domain view.
"""
Base.parentindices(d::SubDomain) = d.inds

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, d::SubDomain{Dim,T}) where {Dim,T}
  name = prettyname(d.domain)
  nelms = length(d.inds)
  print(io, "$nelms view(::$name{$Dim,$T}, ")
  printinds(io, d.inds)
  print(io, ")")
end
