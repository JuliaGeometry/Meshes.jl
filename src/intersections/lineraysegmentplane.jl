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
function intersection(f, g::LineLike{T}, p::Plane{T}) where {T}
  p₀ = origin(p)
  n  = normal(p)

  # get the origin and direction of geometry g
  g₀ = g(0)
  d = g(1) - g(0)
  
  # evaluate denominator
  a = d ⋅ n
  
  # if a is zero, then g is parallel to p
  if isapprox(a, zero(T), atol=atol(T))
    # if the numerator is zero, then g is overlapping
    if isapprox((p₀ - g₀) ⋅ n, zero(T), atol=atol(T))
      return @IT _overlapping(g) g f
    else
      return @IT NoIntersection nothing f
    end
  else
    # calculate the length parameter
    λ = -(n ⋅ (g₀ - p₀)) / a

    return _intersection(f, g, λ)
  end
end