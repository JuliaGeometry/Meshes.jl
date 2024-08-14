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
struct SubDomain{M<:Manifold,C<:CRS,D<:Domain{M,C},I<:AbstractVector{Int}} <: Domain{M,C}
  domain::D
  inds::I
end

# specialize constructor to avoid infinite loops
SubDomain(d::SubDomain, inds::AbstractVector{Int}) = SubDomain(d.domain, d.inds[inds])

"""
    SubGrid{M,CRS,Dim}

A subgrid of geometries in a given manifold `M` with point coordinates specified
in a coordinate reference system `CRS`, which is embedded in `Dim` dimensions.
"""
const SubGrid{M<:Manifold,C<:CRS,Dim} = SubDomain{M,C,<:Grid{M,C,Dim}}

# -----------------
# DOMAIN INTERFACE
# -----------------

element(d::SubDomain, ind::Int) = element(d.domain, d.inds[ind])

nelements(d::SubDomain) = length(d.inds)

# specializations
Base.eltype(d::SubDomain) = eltype(d.domain)

function Base.vcat(d1::SubDomain, d2::SubDomain)
  if d1.domain === d2.domain
    SubDomain(d1.domain, vcat(d1.inds, d2.inds))
  else
    GeometrySet(vcat(collect(d1), collect(d2)))
  end
end

function ==(d1::SubDomain, d2::SubDomain)
  if d1.domain == d2.domain
    d1.inds == d2.inds
  else
    nelements(d1) == nelements(d2) && all(d1[i] == d2[i] for i in 1:nelements(d1))
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

function Base.summary(io::IO, d::SubDomain)
  name = prettyname(d.domain)
  nelm = length(d.inds)
  print(io, "$nelm view(::$name, ")
  printinds(io, d.inds)
  print(io, ")")
end
