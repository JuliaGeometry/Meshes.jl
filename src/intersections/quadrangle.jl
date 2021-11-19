# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(f, segment_or_ray, quadrangle)

Calculate the intersection type of `Segment` or `Ray` and `Quadrangle` and apply function
`f` to it.
"""
function intersecttype(
  f::Function, segment_or_ray::Union{Segment{3,T},Ray{3,T}}, n::Quadrangle{3,T}
) where {T}
  v = vertices(n)
  # Check concave condition
  if (v[1] in Triangle(v[2], v[3], v[4])) || (v[3] in Triangle(v[1], v[2], v[4]))
    intersect1 = segment_or_ray ∩ Triangle(v[1], v[2], v[3])
    intersect2 = segment_or_ray ∩ Triangle(v[1], v[3], v[4])
  else
    intersect1 = segment_or_ray ∩ Triangle(v[1], v[2], v[4])
    intersect2 = segment_or_ray ∩ Triangle(v[2], v[3], v[4])
  end
  # Quadrangle should be in a plane. There is only 1 intersection
  if isnothing(intersect1)
    if isnothing(intersect2)
      return NoIntersection() |> f
    else
      return IntersectingRaySegmentQuadrangle(intersect2) |> f
    end
  else
    return IntersectingRaySegmentQuadrangle(intersect1) |> f
  end
end
