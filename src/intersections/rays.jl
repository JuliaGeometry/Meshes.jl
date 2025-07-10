# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of six types:
#
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at origin of one ray (EdgeTouching -> Point)
# 3. intersect at origin of both rays (CornerTouching -> Point)
# 4. overlap with aligned vectors (PosOverlapping -> Ray)
# 5. overlap with colliding vectors (NegOverlapping -> Segment)
# 6. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, ray₁::Ray, ray₂::Ray)
  a, b = ray₁(0), ray₁(1)
  c, d = ray₂(0), ray₂(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  l₁ = ustrip(norm(b - a))
  l₂ = ustrip(norm(d - c))
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NotIntersecting nothing f #CASE 6
  # collinear
  elseif r == rₐ == 1
    if isnonnegative((b - a) ⋅ (d - c)) # rays aligned in same direction
      if isnonnegative((a - c) ⋅ (b - a)) # origin of ray₁ ∈ ray₂
        return @IT PosOverlapping ray₁ f # CASE 4: ray₁
      else
        return @IT PosOverlapping ray₂ f # CASE 4: ray₂
      end
    else # colliding rays
      if a ∉ ray₂
        return @IT NotIntersecting nothing f # CASE 6
      elseif a == c
        return @IT CornerTouching a f # CASE 3
      else
        return @IT NegOverlapping Segment(a, c) f # CASE 5
      end
    end
    # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(λ₁))
    λ₂ = mayberound(λ₂, zero(λ₂))
    if λ₁ < 0 || λ₂ < 0
      return @IT NotIntersecting nothing f # CASE 6
    elseif λ₁ == 0
      if λ₂ == 0
        return @IT CornerTouching a f # CASE 3
      else
        return @IT EdgeTouching a f # CASE 2: origin of ray₁
      end
    else
      if λ₂ == 0
        return @IT EdgeTouching c f # CASE 2: origin of ray₂
      else
        return @IT Crossing ray₁(λ₁ / l₁) f # CASE 1: equal to ray₂(λ₂/l₂)
      end
    end
  end
end

# The intersection type can be one of four types:
#
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at origin of ray (Touching -> Point)
# 3. overlap of line and ray (Overlapping -> Ray)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, ray::Ray, line::Line)
  a, b = ray(0), ray(1)
  c, d = line(0), line(1)

  # rescaling of point b necessary to gain a parameter λ₁ representing the arc length
  l₁ = ustrip(norm(b - a))
  b₀ = a + 1 / l₁ * (b - a)

  λ₁, _, r, rₐ = intersectparameters(a, b₀, c, d)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NotIntersecting nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT Overlapping ray f # CASE 3
  else # in same plane, not parallel
    λ₁ = mayberound(λ₁, zero(λ₁))
    if λ₁ > 0
      return @IT Crossing ray(λ₁ / l₁) f # CASE 1
    elseif λ₁ == 0
      return @IT Touching a f # CASE 2
    else
      return @IT NotIntersecting nothing f # CASE 4
    end
  end
end

# Williams A, Barrus S, Morley R K, et al., 2005.
# (https://dl.acm.org/doi/abs/10.1145/1198555.1198748)
function intersection(f, ray::Ray, box::Box)
  ℒ = lentype(ray)
  invdir = inv.(ray(1) - ray(0))
  lo, up = to.(extrema(box))
  orig = to(ray(0))

  T = numtype(ℒ)
  tmin = zero(T)
  tmax = typemax(T)

  # check for intersection with slabs along with each axis
  for i in 1:embeddim(ray)
    imin = (lo[i] - orig[i]) * invdir[i]
    imax = (up[i] - orig[i]) * invdir[i]

    # swap variables if necessary
    iinv = invdir[i]
    iinv < zero(iinv) && ((imin, imax) = (imax, imin))

    # the ray is on a face of the box, avoid NaN
    (isnan(imin) || isnan(imax)) && continue

    (tmin > imax || imin > tmax) && return @IT NotIntersecting nothing f

    tmin = max(tmin, imin)
    tmax = min(tmax, imax)
  end

  tmin ≈ tmax && return @IT Touching ray(tmin) f

  return @IT Crossing Segment(ray(tmin), ray(tmax)) f
end

# The intersection type can be one of six types:
#
# 1. Touching - origin of ray intersects middle of triangle
# 2. EdgeTouching - origin of ray intersects edge of triangle
# 3. CornerTouching - origin of ray intersects corner of triangle
# 4. Crossing - middle of ray intersects middle of triangle
# 5. EdgeCrossing - middle of ray intersects edge of triangle
# 6. CornerCrossing - middle of ray intersects corner of triangle
#
# The implementation follows the notation of the non-culling branch of
# Möller, T. & Trumbore, B., 1997 (https://www.tandfonline.com/doi/abs/10.1080/10867651.1997.10487468)
function intersection(f, ray::Ray, tri::Triangle)
  O = ray(0)
  D = ray(1) - ray(0)
  V = vertices(tri)

  E₁ = V[2] - V[1]
  E₂ = V[3] - V[1]
  P = D × E₂

  det = E₁ ⋅ P

  if abs(det) < atol(det)
    # ray is parallel to the plane of the triangle
    return @IT NotIntersecting nothing f
  end

  det⁻¹ = inv(det)

  T = O - V[1]

  # calculate u parameter and test bounds
  u = (T ⋅ P) * det⁻¹
  if u < zero(u) || u > one(u)
    return @IT NotIntersecting nothing f
  end

  Q = T × E₁

  # calculate v parameter and test bounds
  v = (D ⋅ Q) * det⁻¹
  if v < zero(v) || u + v > one(u)
    return @IT NotIntersecting nothing f
  end

  # calculate t parameter and test bounds
  t = (E₂ ⋅ Q) * det⁻¹
  if t < zero(t)
    return @IT NotIntersecting nothing f
  end

  # calculate barycentric coordinates
  w = (u, v, one(u) - u - v)

  if any(isapprox(O), V)
    return @IT CornerTouching ray(t) f
  elseif isapproxzero(t)
    if all(ispositive, w)
      return @IT Touching ray(t) f
    else
      return @IT EdgeTouching ray(t) f
    end
  end

  if count(isapproxzero, w) == 1
    return @IT EdgeCrossing ray(t) f
  elseif count(isapproxone, w) == 1
    return @IT CornerCrossing ray(t) f
  end

  return @IT Crossing ray(t) f
end

intersection(f, ray::Ray, p::Polygon) = intersection(f, GeometrySet([ray]), simplexify(p))
