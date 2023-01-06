# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
Return appropriate type for a geometry overlapping with a `Plane`.
Ideally this would be a macro, but the geometry type isn't known at parse time,
so it is implemented via multiple-dispatch instead
=#
_overlapping(::Line) = OverlappingLinePlane
_overlapping(::Ray) = OverlappingRayPlane
_overlapping(::Segment) = OverlappingSegmentPlane

#=
Intersection of a `Plane` and non-parallel `Line` given a line parameter `λ`.
As a line is infinite, a `CrossingLinePlane` will always be returned.
=#
_intersection(f, l::Line, λ) = @IT CrossingLinePlane l(λ) f

#=
Intersection of a `Plane` and non-parallel `Ray` given ray parameter `λ`.

Return types:
  λ < 0 ⟹ NoIntersection
  λ ≈ 0 ⟹ TouchingRayPlane
  λ > 0 ⟹ CrossingRayPlane
=#
function _intersection(f, r::Ray{3,T}, λ) where {T}
  # if λ is approximately 0, set as so to prevent any domain errors
  if isapprox(λ, zero(T), atol=atol(T))
    return @IT TouchingRayPlane r(zero(T)) f
  end

  # if λ is out of bounds for the ray, then there is no intersection
  if (λ < zero(T))
    return @IT NoIntersection nothing f
  else
    return @IT CrossingRayPlane r(λ) f
  end
end

#=
Intersection of a `Plane` and non-parallel `Segment` given a segment parameter `λ`.

Return types:
  λ < 0 or  λ > 1 ⟹ NoIntersection
  λ ≈ 0 or  λ ≈ 1 ⟹ TouchingSegmentPlane
  λ > 0 and λ < 1 ⟹ CrossingSegmentPlane
=#
function _intersection(f, s::Segment{3,T}, λ) where {T}
  # if λ is approximately 0, set as so to prevent any domain errors
  if isapprox(λ, zero(T), atol=atol(T))
    return @IT TouchingSegmentPlane s(zero(T)) f
  end

  # if λ is approximately 1, set as so to prevent any domain errors
  if isapprox(λ, one(T), atol=atol(T))
    return @IT TouchingSegmentPlane s(one(T)) f
  end

  # if λ is out of bounds for the segment, then there is no intersection
  if (λ < zero(T) || λ > one(T))
    return @IT NoIntersection nothing f
  else
    return @IT CrossingSegmentPlane s(λ) f
  end
end

#=
(https://en.wikipedia.org/wiki/Line-plane_intersection)
=#
const LineLike{T} = Union{Line{3,T}, Ray{3,T}, Segment{3,T}}
function intersection(f, l::LineLike{T}, p::Plane{T}) where {T}
  # origin and direction of line
  l₀ = l(0)
  d  = l(1) - l(0)
  
  # origin and normal of plane
  p₀ = origin(p)
  n  = normal(p)
  
  # auxiliary parameters
  a = (p₀ - l₀) ⋅ n
  b = d ⋅ n
  
  if isapprox(b, zero(T), atol=atol(T))
    if isapprox(a, zero(T), atol=atol(T))
      return @IT _overlapping(l) l f
    else
      return @IT NoIntersection nothing f
    end
  else
    return _intersection(f, l, a / b)
  end
end