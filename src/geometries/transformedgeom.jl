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
  p = t(centroid(g))
  TransformedGeometry{manifold(p),crs(p)}(g, t)
end

# specialize constructor to avoid deep structures
TransformedGeometry(g::TransformedGeometry, t::Transform) = TransformedGeometry(g.geometry, m.transform → t)

# type aliases for convenience
const TransformedPoint{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Point,T}
const TransformedSegment{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Segment,T}
const TransformedRope{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Rope,T}
const TransformedRing{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Ring,T}
const TransformedPolytope{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polytope,T}
const TransformedPolygon{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polygon,T}
const TransformedPolyhedron{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polyhedron,T}

Base.parent(g::TransformedGeometry) = g.geometry

transform(g::TransformedGeometry) = g.transform

# ---------
# GEOMETRY
# ---------

paramdim(g::TransformedGeometry) = paramdim(g.geometry)

centroid(g::TransformedGeometry) = g.transform(centroid(g.geometry))

==(g₁::TransformedGeometry, g₂::TransformedGeometry) = g₁.transform == g₂.transform && g₁.geometry == g₂.geometry

Base.isapprox(g₁::TransformedGeometry, g₂::TransformedGeometry; atol=atol(lentype(g₁)), kwargs...) =
  isapprox(g₁.geometry, g₂.geometry; atol, kwargs...) && g₁.transform == g₂.transform

(g::TransformedGeometry)(uvw...) = g.transform(g.geometry(uvw...))

# ---------
# POLYTOPE
# ---------

vertex(g::TransformedPolytope, ind) = g.transform(vertex(g.geometry, ind))

vertices(g::TransformedPolytope) = map(g.transform, vertices(g.geometry))

nvertices(g::TransformedPolytope) = nvertices(g.geometry)

Base.unique(g::TransformedPolytope) = unique!(deepcopy(g))

Base.unique!(g::TransformedPolytope) = (unique!(g.geometry); g)

# --------
# POLYGON
# --------

rings(p::TransformedPolygon) = map(p.transform, rings(p.geometry))

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, g::TransformedGeometry)
  name = prettyname(g.geometry)
  print(io, "Transformed$name")
end
