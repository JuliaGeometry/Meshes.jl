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

vertex(m::MultiPolytope, ind) = first(Iterators.drop(eachvertex(m), ind - 1))

vertices(m::MultiPolytope) = collect(eachvertex(m))

nvertices(m::MultiPolytope) = sum(nvertices, m.geoms)

eachvertex(m::MultiPolytope) = VertexItr(m)

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

# -----------------
# HELPER FUNCTIONS
# -----------------

struct VertexItr{T}
  el::T
end

_v_iterate(el::Polytope, i) =
  (@inline; (i - 1) % UInt < nvertices(el) % UInt ? (@inbounds vertex(el, i), i + 1) : nothing)

Base.iterate(itr::VertexItr{<:MultiPolytope}, state=(1, 1)) = begin
  ig, ivg = state
  ig > length(itr.el.geoms) && return nothing
  is = _v_iterate(itr.el.geoms[ig], ivg)
  is === nothing && return Base.iterate(itr, (ig + 1, 1))
  v, ivg = is
  return (v, (ig, ivg))
end

Base.IteratorSize(::VertexItr) = Base.HasLength()
Base.IteratorEltype(::VertexItr) = Base.HasEltype()
Base.length(itr::VertexItr{<:MultiPolytope}) = sum(nvertices, itr.el.geoms)
Base.eltype(::VertexItr{<:MultiPolytope{K,M,C}}) where {K,M<:Manifold,C<:CRS} = Point{M,C}
