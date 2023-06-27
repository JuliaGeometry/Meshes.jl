# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const LineLike{T} = Union{Line{3,T},Ray{3,T},Segment{3,T}}

# (https://en.wikipedia.org/wiki/Line-plane_intersection)
function intersection(f, l::LineLike{T}, p::Plane{T}) where {T}
  # origin and direction of line
  l₀ = l(0)
  d = l(1) - l(0)

  # origin and normal of plane
  p₀ = p(0, 0)
  n = normal(p)

  # auxiliary parameters
  a = (p₀ - l₀) ⋅ n
  b = d ⋅ n

  if isapprox(b, zero(T), atol=atol(T))
    if isapprox(a, zero(T), atol=atol(T))
      return @IT Overlapping l f
    else
      return @IT NotIntersecting nothing f
    end
  else
    return _intersection(f, l, a / b)
  end
end

# Intersection of a `Plane` and non-parallel `Line` given a line parameter `λ`.
# As a line is infinite, a `Crossing` will always be returned.
_intersection(f, l::Line, λ) = @IT Crossing l(λ) f

# Intersection of a `Plane` and non-parallel `Ray` given ray parameter `λ`.
# 
# Return types:
#   λ < 0 ⟹ NotIntersecting
#   λ ≈ 0 ⟹ Touching
#   λ > 0 ⟹ Crossing
function _intersection(f, r::Ray{3,T}, λ) where {T}
  # if λ is approximately 0, set as so to prevent any domain errors
  if isapprox(λ, zero(T), atol=atol(T))
    return @IT Touching r(zero(T)) f
  end

  # if λ is out of bounds for the ray, then there is no intersection
  if (λ < zero(T))
    return @IT NotIntersecting nothing f
  else
    return @IT Crossing r(λ) f
  end
end

# Intersection of a `Plane` and non-parallel `Segment` given a segment parameter `λ`.
# 
# Return types:
#   λ < 0 or  λ > 1 ⟹ NotIntersecting
#   λ ≈ 0 or  λ ≈ 1 ⟹ Touching
#   λ > 0 and λ < 1 ⟹ Crossing
function _intersection(f, s::Segment{3,T}, λ) where {T}
  # if λ is approximately 0, set as so to prevent any domain errors
  if isapprox(λ, zero(T), atol=atol(T))
    return @IT Touching s(zero(T)) f
  end

  # if λ is approximately 1, set as so to prevent any domain errors
  if isapprox(λ, one(T), atol=atol(T))
    return @IT Touching s(one(T)) f
  end

  # if λ is out of bounds for the segment, then there is no intersection
  if (λ < zero(T) || λ > one(T))
    return @IT NotIntersecting nothing f
  else
    return @IT Crossing s(λ) f
  end
end
