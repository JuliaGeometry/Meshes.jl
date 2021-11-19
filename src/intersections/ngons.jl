# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(f, segment_or_ray, triangle)

Calculate the intersection type of `segment` or `ray` and `ngon` and apply function `f` to
it. Treat `obj` as a series of triangles.
"""
function intersecttype(
  f::Function, segment_or_ray::Union{Segment{3,T},Ray{3,T}}, n::Ngon{N,3,T}
) where {N,T}
  vs = vertices(n)
  for i in 2:N-1
    intersection = segment_or_ray ∩ Triangle(vs[1], vs[i], vs[i+1])
    if !isnothing(intersection)
      return IntersectingRaySegmentNgon(intersection) |> f
    end
  end
  return NoIntersection() |> f
end
