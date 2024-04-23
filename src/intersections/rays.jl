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
function intersection(f, ray₁::Ray{N}, ray₂::Ray{N}) where {N}
  a, b = ray₁(0), ray₁(1)
  c, d = ray₂(0), ray₂(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  l₁, l₂ = norm(b - a), norm(d - c)
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)
  T = typeof(λ₁)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NotIntersecting nothing f #CASE 6
  # collinear
  elseif r == rₐ == 1
    if (b - a) ⋅ (d - c) ≥ 0 # rays aligned in same direction
      if (a - c) ⋅ (b - a) ≥ 0 # origin of ray₁ ∈ ray₂
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
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(λ₂, zero(T))
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
function intersection(f, ray::Ray{N}, line::Line{N}) where {N}
  a, b = ray(0), ray(1)
  c, d = line(0), line(1)

  # rescaling of point b necessary to gain a parameter λ₁ representing the arc length
  l₁ = norm(b - a)
  b₀ = a + 1 / l₁ * (b - a)

  λ₁, _, r, rₐ = intersectparameters(a, b₀, c, d)
  T = typeof(λ₁)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NotIntersecting nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT Overlapping ray f # CASE 3
  else # in same plane, not parallel
    λ₁ = mayberound(λ₁, zero(T))
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
function intersection(f, ray::Ray{Dim}, box::Box{Dim}) where {Dim}
  invdir = 1 ./ (ray(1) - ray(0))
  lo, up = coordinates.(extrema(box))
  orig = coordinates(ray(0))

  𝒬 = eltype(orig)
  tmin = zero(𝒬)
  tmax = typemax(𝒬)

  # check for intersection with slabs along with each axis
  for i in 1:Dim
    imin = (lo[i] - orig[i]) * invdir[i]
    imax = (up[i] - orig[i]) * invdir[i]

    # swap variables if necessary
    invdir[i] < zero(𝒬) && ((imin, imax) = (imax, imin))

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
# Möller, T. & Trumbore, B., 1997.
# (https://www.tandfonline.com/doi/abs/10.1080/10867651.1997.10487468)
function intersection(f, ray::Ray{3}, tri::Triangle{3})
  vs = vertices(tri)
  o = ray(0)
  d = ray(1) - ray(0)

  e₁ = vs[3] - vs[1]
  e₂ = vs[2] - vs[1]
  p = d × e₂
  det = e₁ ⋅ p

  # keep det > 0, modify T accordingly
  if det > atol(det)
    τ = o - vs[1]
  else
    τ = vs[1] - o
    det = -det
  end

  if det < atol(det)
    # This ray is parallel to the plane of the triangle.
    return @IT NotIntersecting nothing f
  end

  # calculate u parameter and test bounds
  u = τ ⋅ p
  if u < -atol(u) || u > det
    return @IT NotIntersecting nothing f
  end

  q = τ × e₁

  # calculate v parameter and test bounds
  v = d ⋅ q
  if v < -atol(v) || u + v > det
    return @IT NotIntersecting nothing f
  end

  λ = (e₂ ⋅ q) / det

  if λ < -atol(λ)
    return @IT NotIntersecting nothing f
  end

  # assemble barycentric weights
  w = Vec(u, v, det - u - v)

  if any(isapprox.(o, vs))
    return @IT CornerTouching ray(λ) f
  elseif isapprox(λ, zero(λ), atol=atol(λ))
    if all(x -> x > zero(x), w)
      return @IT Touching ray(λ) f
    else
      return @IT EdgeTouching ray(λ) f
    end
  end

  if count(x -> isapprox(x, zero(x), atol=atol(x)), w) == 1
    return @IT EdgeCrossing ray(λ) f
  elseif count(x -> isapprox(x, det, atol=atol(x)), w) == 1
    return @IT CornerCrossing ray(λ) f
  end

  λ = clamp(λ, zero(λ), typemax(λ))

  return @IT Crossing ray(λ) f
end

intersection(f, ray::Ray, p::Polygon) = intersection(f, GeometrySet([ray]), simplexify(p))
