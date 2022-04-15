# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersection(ray, box)

Compute the intersection between a `ray` and a `box`.

## References

* Williams A, Barrus S, Morley R K, et al., 2005. [An efficient and robust ray-box
  intersection algorithm]
  (https://dl.acm.org/doi/abs/10.1145/1198555.1198748)
"""
function intersection(r::Ray{Dim,T}, b::Box{Dim,T}) where {Dim,T}
  invdir = one(T) ./ direction(r)
  lo, up = coordinates.(extrema(b))
  orig = coordinates(origin(r))

  tmin = zero(T)
  tmax = typemax(T)

  # check for intersection with slabs along with each axis
  for i in 1:Dim
    imin = (lo[i] - orig[i]) * invdir[i]
    imax = (up[i] - orig[i]) * invdir[i]

    # swap variables if necessary
    invdir[i] < zero(T) && ((imin, imax) = (imax, imin))

    # the ray is on a face of the box, avoid NaN
    (isnan(imin) || isnan(imax)) && continue

    (tmin > imax || imin > tmax) && return @IT NoIntersection nothing
    
    tmin = max(tmin, imin)
    tmax = min(tmax, imax)
  end

  tmin ≈ tmax && return @IT TouchingRayBox r(tmin)
  
  return @IT CrossingRayBox Segment(r(tmin), r(tmax))
end