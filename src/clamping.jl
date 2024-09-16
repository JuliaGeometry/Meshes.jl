# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    clamp(point, box)

Clamp the coordinates of a [`Point`](@ref) to the edges of a [`Box`](@ref).

For each dimension, coordinates outside of the box are moved to the nearest
edge of the box. The point and box must have an equal number of dimensions.
"""
function Base.clamp(point::Point, box::Box)
  x = cartvalues(point)
  lo = cartvalues(minimum(box))
  hi = cartvalues(maximum(box))
  vals = ntuple(embeddim(point)) do i
    clamp(x[i], lo[i], hi[i])
  end
  Point(withcrs(point, vals))
end

"""
    clamp(pset, box)

Clamp each point in a [`PointSet`](@ref) to the edges of a [`Box`](@ref),
returning a new set of points.
"""
Base.clamp(points::PointSet, box::Box) = PointSet(map(point -> clamp(point, box), points))
