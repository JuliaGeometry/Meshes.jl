# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(f, b1, b2)

Compute the intersection type of two boxes `b1` and `b2`
and apply function `f` to it.

The intersection type can be one of four types:

1. overlap with non-zero measure
2. intersect at one of the boundaries
3. intersect at corner point
4. do not overlap nor intersect
"""
function intersecttype(f::Function, b1::Box{Dim,T}, b2::Box{Dim,T}) where {Dim,T}
  m1, M1 = coordinates.(extrema(b1))
  m2, M2 = coordinates.(extrema(b2))

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  # branch on possible configurations
  if u ≺ v
    return OverlappingBoxes(Box(u, v)) |> f
  elseif u ≻ v
    return NoIntersection() |> f
  elseif isapprox(u, v, atol=atol(T))
    return CornerTouchingBoxes(u) |> f
  else
    return FaceTouchingBoxes(Box(u, v)) |> f
  end
end

"""
    intersecttype(function, ray, box)

Calculate the intersection type of a `ray`` and a `box`` and apply function `f` to it.

## References

* Williams A, Barrus S, Morley R K, et al., 2005. [An efficient and robust ray-box
  intersection algorithm]
  (https://dl.acm.org/doi/abs/10.1145/1198555.1198748)
"""
function intersecttype(f::Function, r::Ray{3,T}, b::Box{3,T}) where {T}
  sign = (x -> x < 0 ? 2 : 1).(r.v)
  invdir = one(T) ./ r.v
  bounds = (b.min, b.max)

  tmin = maximum((bounds[sign[1]] - r.p) .* invdir)
  tmin = max(tmin, zero(T))
  tmax = minimum((bounds[3 - sign[1]] - r.p) .* invdir)

  if tmin < tmax
    return CrossingRayBox(Segment(r(tmin), r(tmax))) |> f
  elseif tmin ≈ tmax
    return TouchingRayBox(r(tmin)) |> f
  else
    return NoIntersection() |> f
  end
end
