# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Collection(items)

A collection of `items` (points or geometries) seen as a single domain.
"""
struct Collection{Dim,T,I<:PointOrGeometry{Dim,T}} <: Domain{Dim,T}
  items::Vector{I}
end

# -----------------
# DOMAIN INTERFACE
# -----------------

element(c::Collection, ind::Int) = c.items[ind]

nelements(c::Collection) = length(c.items)

# ------------------------
# SPECIAL CASE: POINT SET
# ------------------------

const PointSet{Dim,T} = Collection{Dim,T,Point{Dim,T}}

"""
    PointSet(points)

A set of `points` (a.k.a. point cloud) seen as a single domain.

## Example

All point sets below are the same and contain two points in R³:

```julia
julia> PointSet([Point(1,2,3), Point(4,5,6)])
julia> PointSet(Point(1,2,3), Point(4,5,6))
julia> PointSet([(1,2,3), (4,5,6)])
julia> PointSet((1,2,3), (4,5,6))
julia> PointSet([[1,2,3], [4,5,6]])
julia> PointSet([1,2,3], [4,5,6])
julia> PointSet([1 4; 2 5; 3 6])
```
"""
PointSet(points::AbstractVector{P}) where {P<:Point} = PointSet{embeddim(P),coordtype(P)}(points)
PointSet(points::Vararg{P}) where {P<:Point} = PointSet(collect(points))
PointSet(coords::AbstractVector{TP}) where {TP<:Tuple} = PointSet(Point.(coords))
PointSet(coords::Vararg{TP}) where {TP<:Tuple} = PointSet(collect(coords))
PointSet(coords::AbstractVector{V}) where {V<:AbstractVector} = PointSet(Point.(coords))
PointSet(coords::Vararg{V}) where {V<:AbstractVector} = PointSet(collect(coords))
PointSet(coords::AbstractMatrix) = PointSet(Point.(eachcol(coords)))

centroid(pset::PointSet, ind::Int) = pset[ind]

function Base.show(io::IO, pset::PointSet{Dim,T}) where {Dim,T}
  nelm = nelements(pset)
  print(io, "$nelm PointSet{$Dim,$T}")
end

# ---------------------------
# SPECIAL CASE: GEOMETRY SET
# ---------------------------

const GeometrySet{Dim,T,G<:Geometry{Dim,T}} = Collection{Dim,T,G}

"""
    GeometrySet(geometries)

A set of `geometries` seen as a single domain.
"""
GeometrySet(geometries::AbstractVector{G}) where {G<:Geometry} =
  GeometrySet{embeddim(G),coordtype(G),G}(geometries)

function Base.show(io::IO, gset::GeometrySet{Dim,T}) where {Dim,T}
  nelm = nelements(gset)
  print(io, "$nelm GeometrySet{$Dim,$T}")
end
