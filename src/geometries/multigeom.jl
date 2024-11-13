# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Multi(geoms)

A collection of geometries `geoms` seen as a single [`Geometry`](@ref).

In geographic information systems (GIS) it is common to represent
multiple polygons as a single entity (e.g. country with islands).

### Notes

- Type aliases are [`MultiPoint`](@ref), [`MultiSegment`](@ref),
  [`MultiRope`](@ref), [`MultiRing`](@ref), [`MultiPolygon`](@ref).
"""
struct Multi{M<:Manifold,C<:CRS,G<:Geometry{M,C}} <: Geometry{M,C}
  geoms::Vector{G}
end

# constructor with iterator of geometries
Multi(geoms) = Multi(collect(geoms))

# type aliases for convenience
const MultiPoint{M<:Manifold,C<:CRS} = Multi{M,C,<:Point{M,C}}
const MultiSegment{M<:Manifold,C<:CRS} = Multi{M,C,<:Segment{M,C}}
const MultiRope{M<:Manifold,C<:CRS} = Multi{M,C,<:Rope{M,C}}
const MultiRing{M<:Manifold,C<:CRS} = Multi{M,C,<:Ring{M,C}}
const MultiPolygon{M<:Manifold,C<:CRS} = Multi{M,C,<:Polygon{M,C}}
const MultiPolyhedron{M<:Manifold,C<:CRS} = Multi{M,C,<:Polyhedron{M,C}}
const MultiPolytope{K,M<:Manifold,C<:CRS} = Multi{M,C,<:Polytope{K,M,C}}

Base.parent(m::Multi) = m.geoms

# ---------
# GEOMETRY
# ---------

paramdim(m::Multi) = maximum(paramdim, m.geoms)

==(m₁::Multi, m₂::Multi) = m₁.geoms == m₂.geoms

Base.isapprox(m₁::Multi, m₂::Multi; atol=atol(lentype(m₁)), kwargs...) =
  length(m₁.geoms) == length(m₂.geoms) && all(isapprox(g₁, g₂; atol, kwargs...) for (g₁, g₂) in zip(m₁.geoms, m₂.geoms))

# ---------
# POLYTOPE
# ---------

vertex(m::MultiPolytope, ind) = first(Base.Iterators.drop(eachvertex(m), ind - 1)) # nth(itr, n)

vertices(m::MultiPolytope) = collect(eachvertex(m))

nvertices(m::MultiPolytope) = sum(nvertices, m.geoms)

eachvertex(m::MultiPolytope) = (v for g in m.geoms for v in vertices(g))

Base.unique(m::MultiPolytope) = unique!(deepcopy(m))

function Base.unique!(m::MultiPolytope)
  foreach(unique!, m.geoms)
  m
end

# --------
# POLYGON
# --------

rings(m::MultiPolygon) = [ring for poly in m.geoms for ring in rings(poly)]

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, m::Multi)
  name = prettyname(eltype(m.geoms))
  print(io, "Multi$name")
end

function Base.show(io::IO, m::Multi)
  print(io, "Multi(")
  geoms = prettyname.(m.geoms)
  counts = ("$(count(==(g), geoms))×$g" for g in unique(geoms))
  join(io, counts, ", ")
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::Multi)
  summary(io, m)
  println(io)
  printelms(io, m.geoms)
end
