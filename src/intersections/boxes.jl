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
function intersecttype(f::Function, r::Ray{Dim,T}, b::Box{Dim,T}) where {Dim,T}
  invdir = one(T) ./ direction(r)
  lo, up = coordinates.(extrema(b))
  orig = coordinates(origin(r))

  tmin = zero(T) / oneunit(T)
  tmax = typemax(T) / oneunit(T)

  # check for intersection with slabs along with each axis
  for i in 1:Dim
    imin = (lo[i] - orig[i]) * invdir[i]
    imax = (up[i] - orig[i]) * invdir[i]

    # swap variables if necessary
    invdir[i] < zero(eltype(invdir)) && ((imin, imax) = (imax, imin))

    # the ray is on a face of the box, avoid NaN
    (isnan(imin) || isnan(imax)) && continue

    (tmin > imax || imin > tmax) && return NoIntersection() |> f
    
    tmin = max(tmin, imin)
    tmax = min(tmax, imax)
  end

  tmin ≈ tmax && return TouchingRayBox(r(tmin)) |> f
  
  return CrossingRayBox(Segment(r(tmin), r(tmax))) |> f
end
