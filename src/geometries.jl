# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{Dim,T}

A geometry embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Geometry{Dim,T} end

"""
    embeddim(geometry)

Return the number of dimensions of the space where the `geometry` is embedded.
"""
embeddim(::Type{<:Geometry{Dim,T}}) where {Dim,T} = Dim
embeddim(g::Geometry) = embeddim(typeof(g))

"""
    paramdim(geometry)

Return the number of parametric dimensions of the `geometry`. For example, a
sphere embedded in 3D has 2 parametric dimension (polar and azimuthal angles).
"""
paramdim(g::Geometry) = paramdim(typeof(g))

"""
    coordtype(geometry)

Return the machine type of each coordinate used to describe the `geometry`.
"""
coordtype(::Type{<:Geometry{Dim,T}}) where {Dim,T} = T
coordtype(g::Geometry) = coordtype(typeof(g))

"""
    centroid(geometry)

Return the centroid of the `geometry`.
"""
centroid(g::Geometry) = center(g)

"""
    measure(geometry)

Return the measure of the `geometry`, i.e. the length, area, or volume.
"""
function measure end

"""
    boundary(geometry)

Return the boundary of the `geometry`.
"""
function boundary end

"""
    g₁ ⊆ g₂

Tells whether or not the geometry `g₁` is a subset of geometry `g₂`.
"""
Base.issubset(g₁::Geometry, g₂::Geometry) = all(p ∈ g₂ for p in vertices(g₁))

"""
    p ∈ g

Tells whether or not the point `p` is in the geometry `g`.
"""
Base.in(p::Point, g::Geometry)

"""
    isconvex(geometry)

Tells whether or not the `geometry` is convex.
"""
isconvex(g::Geometry) = isconvex(typeof(g))

"""
    issimplex(geometry)

Tells whether or not the `geometry` is simplex.
"""
issimplex(::Type{<:Geometry}) = false
issimplex(g::Geometry) = issimplex(typeof(g))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("primitives.jl")
include("polytopes.jl")

"""
    PointOrGeometry{Dim,T}

A union type that can represent either a point or a geometry.
"""
const PointOrGeometry{Dim,T} = Union{Point{Dim,T},Geometry{Dim,T}}

"""
    Multi(entities)

A collection of points or geometries seen as a single entity.

In geographic information systems (GIS) it is common to represent
multiple polygons as a single entity. In this case the polygons
refer to the same object in the real world (e.g. country with islands).
"""
struct Multi{Dim,T,I<:PointOrGeometry{Dim,T}} <: Geometry{Dim,T}
  items::Vector{I}
end

Base.getindex(multi::Multi, ind) = getindex(multi.items, ind)
Base.length(multi::Multi) = length(multi.items)
Base.eltype(multi::Multi) = eltype(multi.items)
Base.firstindex(multi::Multi) = firstindex(multi.items)
Base.lastindex(multi::Multi) = lastindex(multi.items)
Base.iterate(multi::Multi, state=1) =
  state > length(multi) ? nothing : (multi[state], state+1)

paramdim(multi::Multi) = maximum(paramdim, multi.items)

==(multi₁::Multi, multi₂::Multi) =
  length(multi₁) == length(multi₂) &&
  all(g -> g[1] == g[2], zip(multi₁, multi₂))

vertices(multi::Multi) =
  [vertex for geom in multi for vertex in vertices(geom)]

function centroid(multi::Multi)
  cs = coordinates.(centroid.(multi.items))
  Point(sum(cs) / length(cs))
end

measure(multi::Multi) = sum(measure, multi.items)

area(multi::Multi{Dim,T,<:Polygon}) where{Dim,T} = measure(multi)

boundary(multi::Multi) = Multi([boundary(geom) for geom in multi])

chains(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} =
  [chain for geom in multi for chain in chains(geom)]

Base.in(point::Point, multi::Multi) = any(geom -> point ∈ geom, multi)

function Base.show(io::IO, multi::Multi{Dim,T}) where {Dim,T}
  n = length(multi.items)
  G = eltype(multi.items)
  print(io, "$n Multi-$(nameof(G)){$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", multi::Multi)
  println(io, multi)
  print(io, io_lines(multi, "  "))
end
