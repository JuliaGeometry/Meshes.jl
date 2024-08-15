# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TransformedGeometry(geometry, transform)

Lazy representation of a geometric `transform` applied to a `geometry`.
"""
struct TransformedGeometry{M<:Manifold,C<:CRS,G<:Geometry,T<:Transform} <: Geometry{M,C}
  geometry::G
  transform::T

  function TransformedGeometry{M,C}(geometry::G, transform::T) where {M<:Manifold,C<:CRS,G<:Geometry,T<:Transform}
    new{M,C,G,T}(geometry, transform)
  end
end

function TransformedGeometry(g::Geometry, t::Transform)
  p = t(_point(g))
  TransformedGeometry{manifold(p),crs(p)}(g, t)
end

_point(g) = isparametrized(g) ? g(ntuple(i -> zero(numtype(lentype(g))), paramdim(g))...) : centroid(g)

# specialize constructor to avoid deep structures
TransformedGeometry(g::TransformedGeometry, t::Transform) = TransformedGeometry(g.geometry, g.transform → t)

# type aliases for convenience
const TransformedPoint{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Point,T}
const TransformedSegment{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Segment,T}
const TransformedRope{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Rope,T}
const TransformedRing{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Ring,T}
const TransformedPolygon{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polygon,T}
const TransformedPolyhedron{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polyhedron,T}
const TransformedPolytope{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polytope,T}

Base.parent(g::TransformedGeometry) = g.geometry

transform(g::TransformedGeometry) = g.transform

# ---------
# GEOMETRY
# ---------

paramdim(g::TransformedGeometry) = paramdim(g.geometry)

==(g₁::TransformedGeometry, g₂::TransformedGeometry) = pointify(g₁) == pointify(g₂)

==(g₁::TransformedGeometry, g₂::Geometry) = pointify(g₁) == pointify(g₂)

==(g₁::Geometry, g₂::TransformedGeometry) = g₂ == g₁

Base.isapprox(g₁::TransformedGeometry, g₂::TransformedGeometry; atol=atol(lentype(g₁)), kwargs...) =
  all(isapprox(p₁, p₂; atol, kwargs...) for (p₁, p₂) in zip(pointify(g₁), pointify(g₂)))

Base.isapprox(g₁::TransformedGeometry, g₂::Geometry; atol=atol(lentype(g₁)), kwargs...) =
  all(isapprox(p₁, p₂; atol, kwargs...) for (p₁, p₂) in zip(pointify(g₁), pointify(g₂)))

Base.isapprox(g₁::Geometry, g₂::TransformedGeometry; atol=atol(lentype(g₁)), kwargs...) =
  isapprox(g₂, g₁; atol, kwargs...)

(g::TransformedGeometry)(uvw...) = g.transform(g.geometry(uvw...))

# ---------
# POLYTOPE
# ---------

vertex(p::TransformedPolytope, ind) = p.transform(vertex(p.geometry, ind))

vertices(p::TransformedPolytope) = map(p.transform, vertices(p.geometry))

nvertices(p::TransformedPolytope) = nvertices(p.geometry)

Base.unique(p::TransformedPolytope) = unique!(deepcopy(p))

Base.unique!(p::TransformedPolytope) = (unique!(p.geometry); p)

# --------
# POLYGON
# --------

rings(p::TransformedPolygon) = map(p.transform, rings(p.geometry))

# -----------
# IO METHODS
# -----------

prettyname(g::TransformedGeometry) = "Transformed$(prettyname(g.geometry))"
