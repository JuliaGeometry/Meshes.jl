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
measure(g::Geometry) = sum(measure, discretize(g))

"""
    perimeter(geometry)

Return the perimeter of the `geometry`, i.e. the measure of its boundary.
"""
perimeter(g::Geometry) = measure(boundary(g))

"""
    extrema(geometry)

Return the top left and bottom right corners of the
bounding box of the `geometry`.
"""
Base.extrema(g::Geometry) = extrema(boundingbox(g))

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

"""
    isperiodic(geometry)

Tells whether or not the `geometry` is periodic
along each parametric dimension.
"""
isperiodic(g::Geometry) = isperiodic(typeof(g))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("polytopes.jl")
include("primitives.jl")

"""
    PointOrGeometry{Dim,T}

A union type that can represent either a point or a geometry.
"""
const PointOrGeometry{Dim,T} = Union{Point{Dim,T},Geometry{Dim,T}}

"""
    Multi(items)

A collection of points or geometries seen as a single item.

In geographic information systems (GIS) it is common to represent
multiple polygons as a single item. In this case the polygons refer
to the same object in the real world (e.g. country with islands).
"""
struct Multi{Dim,T,I<:PointOrGeometry{Dim,T}} <: Geometry{Dim,T}
  items::Vector{I}
end

# constructor with iterator of items
Multi(items) = Multi(collect(items))

paramdim(multi::Multi) = maximum(paramdim, multi.items)

vertex(multi::Multi, ind) = vertices(multi)[ind]

vertices(multi::Multi) = [vertex for geom in multi.items for vertex in vertices(geom)]

nvertices(multi::Multi) = sum(nvertices, multi.items)

function centroid(multi::Multi)
  cs = coordinates.(centroid.(multi.items))
  Point(sum(cs) / length(cs))
end

measure(multi::Multi) = sum(measure, multi.items)

Base.length(multi::Multi{Dim,T,<:Polytope{1}}) where {Dim,T} = measure(multi)
area(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} = measure(multi)
volume(multi::Multi{Dim,T,<:Polyhedron}) where {Dim,T} = measure(multi)

function boundary(multi::Multi)
  bounds = [boundary(geom) for geom in multi.items]
  valid = filter(!isnothing, bounds)
  isempty(valid) ? nothing : reduce(merge, valid)
end

rings(multi::Multi{Dim,T,<:Polygon}) where {Dim,T} = [ring for poly in multi.items for ring in rings(poly)]

Base.collect(multi::Multi) = multi.items

Base.in(point::Point, multi::Multi) = any(geom -> point ∈ geom, multi.items)

==(multi₁::Multi, multi₂::Multi) =
  length(multi₁.items) == length(multi₂.items) && all(g -> g[1] == g[2], zip(multi₁.items, multi₂.items))

function Base.show(io::IO, multi::Multi{Dim,T}) where {Dim,T}
  name = prettyname(eltype(multi.items))
  print(io, "Multi$name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", multi::Multi)
  println(io, multi)
  print(io, io_lines(multi, "  "))
end
