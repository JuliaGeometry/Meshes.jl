# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometrySet(geometries)

A set of `geometries` representing a [`Domain`](@ref).

## Examples

Set containing two balls centered at `(0.0, 0.0)` and `(1.0, 1.0)`:

```julia
julia> GeometrySet([Ball((0.0, 0.0)), Ball((1.0, 1.0))])
```
"""
struct GeometrySet{Dim,T,G<:Geometry{Dim,T}} <: Domain{Dim,T}
  geoms::Vector{G}
end

# constructor with iterator of geometries
GeometrySet(geoms) = GeometrySet(collect(geoms))

element(gset::GeometrySet, ind::Int) = gset.geoms[ind]

nelements(gset::GeometrySet) = length(gset.geoms)

Base.parent(gset::GeometrySet) = gset.geoms

# ------------------------
# SPECIAL CASE: POINT SET
# ------------------------

const PointSet{Dim,T} = GeometrySet{Dim,T,Point{Dim,T}}

"""
    PointSet(points)

A set of `points` (a.k.a. point cloud) representing a [`Domain`](@ref).

## Examples

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
PointSet(coords::AbstractMatrix) = PointSet(Tuple.(eachcol(coords)))

# constructor with iterator of points
PointSet(points) = PointSet(collect(points))

centroid(pset::PointSet, ind::Int) = pset[ind]

centroid(pset::PointSet) = Point(sum(coordinates, pset) / nelements(pset))
