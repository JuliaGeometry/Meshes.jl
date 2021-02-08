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
struct PointSet{Dim,T}
  points::Vector{Point{Dim,T}} 
end

PointSet(points::Vararg{P}) where {P<:Point} = PointSet(collect(points))
PointSet(points::AbstractVector{TP}) where {TP<:Tuple} = PointSet(Point.(points))
PointSet(points::Vararg{TP}) where {TP<:Tuple} = PointSet(collect(points))

==(pset1::PointSet, pset2::PointSet) = pset1.points == pset2.points

"""
    embeddim(pointset)

Return the number of dimensions of the space where the `pointset` is embedded.
"""
embeddim(::Type{<:PointSet{Dim,T}}) where {Dim,T} = Dim
embeddim(pset::PointSet) = embeddim(typeof(pset))

"""
    coordtype(pointset)

Return the machine type of each coordinate used to describe the `pointset`.
"""
coordtype(::Type{<:PointSet{Dim,T}}) where {Dim,T} = T
coordtype(pset::PointSet) = coordtype(typeof(pset))

"""
    coordinates!(buff, pointset, ind)

Compute the coordinates `buff` of the `ind`-th point in the `pointset` in place.
"""
function coordinates!(buff, pset::PointSet, ind::Int)
  buff .= coordinates(pset.points[ind])
end

# ------------
# IO methods
# ------------
function Base.show(io::IO, pset::PointSet{Dim,T}) where {Dim,T}
  npts = length(pset.points)
  print(io, "$npts PointSet{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", pset::PointSet{Dim,T}) where {Dim,T}
  println(io, pset)
  lines = ["  └─$p" for p in pset.points]
  lines = length(lines) > 11 ? [lines[begin:5]; ["  ⋮"]; lines[end-4:end]] : lines
  print(io, join(lines, "\n"))
end
