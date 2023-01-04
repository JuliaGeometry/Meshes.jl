# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Return appropriate type for a geometry overlapping with a `Plane`

### Notes

- Ideally this would be a macro, but the geometry type isn't known at parse time, so it is implemented with multiple-dispatch instead
"""
function planeintersectionoverlapping(::Line)
  return OverlappingLinePlane
end

function planeintersectionoverlapping(::Ray)
  return OverlappingRayPlane
end

function planeintersectionoverlapping(::Segment)
  return OverlappingSegmentPlane
end

"""
    checkparameterplaneintersection(f, l::Line, λ)

Return the intersection for a `Plane` and non-parallel `Line`, `l` given a line parameter `λ`. As a line is infinite, a `CrossingLinePlane` will always be returned.
"""
function checkparameterplaneintersection(f, l::Line, λ)
  return @IT CrossingLinePlane l(λ) f
end

"""
    checkparameterplaneintersection(f, r::Ray, λ)

Return the intersection for a `Plane` and non-parallel `Ray`, `r` given a ray parameter `λ`.

Return type:
  λ < 0 ⟹ NoIntersection
  λ ≈ 0 ⟹ TouchingRayPlane
  λ > 0 ⟹ CrossingRayPlane
"""
function checkparameterplaneintersection(f, r::Ray{3,T}, λ) where {T}
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

"""
    checkparameterplaneintersection(f, s::Segment, λ)

Return the intersection for a `Plane` and non-parallel `Segment`, `s` given a segment parameter `λ`.

Return type:
  λ < 0 or  λ > 1 ⟹ NoIntersection
  λ ≈ 0 or  λ ≈ 1 ⟹ TouchingSegmentPlane
  λ > 0 and λ < 1 ⟹ CrossingSegmentPlane
"""
function checkparameterplaneintersection(f, s::Segment{3,T}, λ) where {T}
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
function Meshes.intersection(f, g::G, p::Plane{T}) where {T, G<:Union{Line{3,T}, Ray{3,T}, Segment{3,T}}}
  p₀ = coordinates(origin(p))
  n  = normal(p)

  # Get the origin and direction of geometry g
  g₀ = g(0)
  gdir = g(1) - g(0)
  
  # calculate components
  ln = gdir ⋅ n
  
  # if ln is zero, then g is parallel to the plane
  if isapprox(ln, zero(T), atol=atol(T))
    # if the numerator is zero, then g is overlapping
    if isapprox(coordinates(p₀ - g₀) ⋅ n, zero(T), atol=atol(T))
      return @IT planeintersectionoverlapping(g) g f
    else
      return @IT NoIntersection nothing f
    end
  else
    # calculate the length parameter
    λ = -(n ⋅ coordinates(g₀ - p₀)) / ln

    return checkparameterplaneintersection(f, g, λ)
  end
end