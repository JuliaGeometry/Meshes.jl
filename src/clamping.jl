"""
    clamp(point, box)

Clamp the coordinates of a [`Point`](@ref) to the edges of a [`Box`](@ref). For each dimension, coordinates outside of the box are moved to the nearest edge of the box. The point and box must have an equal number of dimensions."""
function Base.clamp(point::Point{Dim,T}, box::Box{Dim,T})::Point{Dim,T} where {Dim,T}
  x = coordinates(point)
  lo = coordinates(minimum(box))
  hi = coordinates(maximum(box))
  ntuple(Dim) do i
    clamp(x[i], lo[i], hi[i])
  end |> Point
end

"""
    clamp(pointset, box)

Clamp each point in a [`PointSet`](@ref) to the edges of a [`Box`](@ref), returning a new set of points."""
function Base.clamp(points::PointSet{Dim,T}, box::Box{Dim,T})::PointSet{Dim,T} where {Dim,T}
  PointSet(map(point -> clamp(point, box), points))
end
