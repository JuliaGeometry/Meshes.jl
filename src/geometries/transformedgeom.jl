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

Base.parent(g::TransformedGeometry) = g.geometry

transform(g::TransformedGeometry) = g.transform

paramdim(g::TransformedGeometry) = paramdim(g.geometry)

centroid(g::TransformedGeometry) = g.transform(centroid(g.geometry))

vertex(g::TransformedGeometry, ind) = vertices(g)[ind]

vertices(g::TransformedGeometry) = map(g.transform, vertices(g.geometry))

nvertices(g::TransformedGeometry) = nvertices(g.geometry)

Base.unique(g::TransformedGeometry) = unique!(deepcopy(g))

function Base.unique!(g::TransformedGeometry)
  unique!(g.geometry)
  g
end

==(g₁::TransformedGeometry, g₂::TransformedGeometry) = g₁.transform == g₂.transform && g₁.geometry == g₂.geometry

Base.isapprox(g₁::TransformedGeometry, g₂::TransformedGeometry; atol=atol(lentype(g₁)), kwargs...) =
  isapprox(g₁.geometry, g₂.geometry; atol, kwargs...) && g₁.transform == g₂.transform

# alias to improve readability in IO methods
const TransformedPolygon{M<:Manifold,C<:CRS,T<:Transform} = TransformedGeometry{M,C,<:Polygon,T}

rings(p::TransformedPolygon) = map(p.transform, rings(p.geometry))

function Base.summary(io::IO, g::TransformedGeometry)
  name = prettyname(g.geometry)
  print(io, "Transformed$name")
end
