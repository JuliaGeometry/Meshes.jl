# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PointSet(points)

A set of `points` (a.k.a. point cloud) seen as a single entity.

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
struct PointSet{Dim,T} <: Domain{Dim,T}
  points::Vector{Point{Dim,T}} 
end

PointSet(points::Vararg{P}) where {P<:Point} = PointSet(collect(points))
PointSet(coords::AbstractVector{TP}) where {TP<:Tuple} = PointSet(Point.(coords))
PointSet(coords::Vararg{TP}) where {TP<:Tuple} = PointSet(collect(coords))
PointSet(coords::AbstractVector{V}) where {V<:AbstractVector} = PointSet(Point.(coords))
PointSet(coords::Vararg{V}) where {V<:AbstractVector} = PointSet(collect(coords))
PointSet(coords::AbstractMatrix) = PointSet(Point.(eachcol(coords)))

==(pset1::PointSet, pset2::PointSet) = pset1.points == pset2.points

# -----------------
# DOMAIN INTERFACE
# -----------------

Base.getindex(pset::PointSet, ind::Int) = getindex(pset.points, ind)

nelements(pset::PointSet) = length(pset.points)

function coordinates!(buff, pset::PointSet, ind::Int)
  buff .= coordinates(pset.points[ind])
end
