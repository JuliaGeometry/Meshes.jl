# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PointSet(points)

A set of `points` (a.k.a. point cloud) seen as a single entity.

## Examples

Create a 2D point set with 100 points:

```julia
julia> PointSet(rand(Point2, 100))
```
"""
struct PointSet{Dim,T} <: Discretization{Dim,T}
  points::Vector{Point{Dim,T}} 
end

PointSet(points::Vararg{P}) where {P<:Point} = PointSet(collect(points))
PointSet(points::AbstractVector{TP}) where {TP<:Tuple} = PointSet(Point.(points))
PointSet(points::Vararg{TP}) where {TP<:Tuple} = PointSet(collect(points))

==(pset1::PointSet, pset2::PointSet) = pset1.points == pset2.points

Base.getindex(pset::PointSet, ind::Int) = getindex(pset.points, ind)

nelements(pset::PointSet) = length(pset.points)

function coordinates!(buff, pset::PointSet, ind::Int)
  buff .= coordinates(pset.points[ind])
end

# -----------
# IO methods
# -----------

function Base.show(io::IO, pset::PointSet{Dim,T}) where {Dim,T}
  npts = length(pset.points)
  print(io, "$npts PointSet{$Dim,$T}")
end
