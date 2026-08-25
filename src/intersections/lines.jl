# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of three types:
# 1. intersect at one point (Crossing -> Point)
# 2. overlap at more than one point (Overlapping -> Line)
# 3. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, line₁::Line, line₂::Line)
  a, b = line₁(0), line₁(1)
  c, d = line₂(0), line₂(1)

  λ₁, _, r, rₐ = intersectparameters(a, b, c, d)

  if r == rₐ == 2
    return @IT Crossing (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT Overlapping line₁ f
  else
    return @IT NotIntersecting nothing f
  end
end

# The intersection type can be one of three types:
# 1. intersect at one point (Touching -> Point)
# 2. intersect at two points (Crossing -> Segment)
# 3. do not intersect (NotIntersecting -> Nothing)
#
# Williams A, Barrus S, Morley R K, et al., 2005.
# (https://dl.acm.org/doi/abs/10.1145/1198555.1198748)
function intersection(f, line::Line, box::Box)
  ℒ = lentype(line)
  invdir = inv.(line(1) - line(0))
  lo, up = to.(extrema(box))
  orig = to(line(0))

  T = numtype(ℒ)
  tmin = typemin(T)
  tmax = typemax(T)

  # check for intersection with slabs along with each axis
  for i in 1:embeddim(line)
    imin = (lo[i] - orig[i]) * invdir[i]
    imax = (up[i] - orig[i]) * invdir[i]

    # swap variables if necessary
    iinv = invdir[i]
    iinv < zero(iinv) && ((imin, imax) = (imax, imin))

    # the line is on a face of the box, avoid NaN
    (isnan(imin) || isnan(imax)) && continue

    (tmin > imax || imin > tmax) && return @IT NotIntersecting nothing f

    tmin = max(tmin, imin)
    tmax = min(tmax, imax)
  end

  tmin ≈ tmax && return @IT Touching line(tmin) f

  return @IT Crossing Segment(line(tmin), line(tmax)) f
end
