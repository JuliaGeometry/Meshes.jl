# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const LineLike{T} = Union{Line{3,T},Ray{3,T},Segment{3,T}}

# (https://en.wikipedia.org/wiki/Line-plane_intersection)
function intersection(f, line::LineLike{T}, plane::Plane{T}) where {T}
  # auxiliary parameters
  d = line(1) - line(0)
  n = normal(plane)
  a = (plane(0, 0) - line(0)) ⋅ n
  b = d ⋅ n
  if isapprox(b, zero(T), atol=atol(T))
    if isapprox(a, zero(T), atol=atol(T))
      return @IT Overlapping line f
    else
      return @IT NotIntersecting nothing f
    end
  else
    return _intersection(f, line, a / b)
  end
end

# (https://en.wikipedia.org/wiki/Plane-plane_intersection)
function intersection(f, plane1::Plane{T}, plane2::Plane{T}) where {T}
  n1 = normal(plane1)
  n2 = normal(plane2)
  n1n2 = n1 ⋅ n2

  if isapprox(n1n2, one(T), atol=atol(T))
    # planes are parallel and do not intersect
    return @IT NotIntersecting nothing f
  else
    d = n1 × n2
    h1 = n1 ⋅ plane1.p.coords
    h2 = n2 ⋅ plane2.p.coords
    c1 = (h1 - h2 * n1n2) / (1 - n1n2^2)
    c2 = (h2 - h1 * n1n2) / (1 - n1n2^2)
    p1 = (c1 * n1) + (c2 * n2)
    p2 = p1 + d
    return @IT Intersecting Line(p1,p2) f
end

# Intersection of a `Plane` and non-parallel `Segment` given a segment parameter `λ`.
# 
# Return types:
#   λ < 0 or  λ > 1 ⟹ NotIntersecting
#   λ ≈ 0 or  λ ≈ 1 ⟹ Touching
#   λ > 0 and λ < 1 ⟹ Crossing
function _intersection(f, seg::Segment{3,T}, λ) where {T}
  # if λ is approximately 0, set as so to prevent any domain errors
  if isapprox(λ, zero(T), atol=atol(T))
    return @IT Touching seg(0) f
  end

  # if λ is approximately 1, set as so to prevent any domain errors
  if isapprox(λ, one(T), atol=atol(T))
    return @IT Touching seg(1) f
  end

  # if λ is out of bounds for the segment, then there is no intersection
  if (λ < zero(T) || λ > one(T))
    return @IT NotIntersecting nothing f
  else
    return @IT Crossing seg(λ) f
  end
end

# Intersection of a `Plane` and non-parallel `Ray` given ray parameter `λ`.
# 
# Return types:
#   λ < 0 ⟹ NotIntersecting
#   λ ≈ 0 ⟹ Touching
#   λ > 0 ⟹ Crossing
function _intersection(f, ray::Ray{3,T}, λ) where {T}
  # if λ is approximately 0, set as so to prevent any domain errors
  if isapprox(λ, zero(T), atol=atol(T))
    return @IT Touching ray(0) f
  end

  # if λ is out of bounds for the ray, then there is no intersection
  if (λ < zero(T))
    return @IT NotIntersecting nothing f
  else
    return @IT Crossing ray(λ) f
  end
end

# Intersection of a `Plane` and non-parallel `Line` given a line parameter `λ`.
# As a line is infinite, a `Crossing` will always be returned.
_intersection(f, line::Line, λ) = @IT Crossing line(λ) f
