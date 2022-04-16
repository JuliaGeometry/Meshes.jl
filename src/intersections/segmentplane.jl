# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
(https://en.wikipedia.org/wiki/Line-plane_intersection)
=#
function intersection(f, s::Segment{3,T}, p::Plane{T}) where {T}
  sᵥ = coordinates.(vertices(s))
  pₒ = coordinates(origin(p))
  n  = normal(p)
  
  # calculate components
  ln = (sᵥ[2] - sᵥ[1]) ⋅ n

  # if ln is zero, the segment is parallel to the plane
  if isapprox(ln, zero(T), atol=atol(T))
    # if the numerator is zero, the segment is coincident
    if isapprox((pₒ - sᵥ[1]) ⋅ n, zero(T), atol=atol(T))
      return @IT OverlappingSegmentPlane s f
    else
      return @IT NoIntersection nothing f
    end
  else
    # calculate the segment parameter
    λ = ((pₒ - sᵥ[1]) ⋅ n) / ln

    # if λ is approximately 0 or 1, set as so to prevent any domain errors
    if isapprox(λ, zero(T), atol=atol(T))
      return @IT TouchingSegmentPlane s(zero(T)) f
    elseif isapprox(λ, one(T), atol=atol(T))
      return @IT TouchingSegmentPlane s(one(T)) f
    end

    # if λ is out of bounds for the segment, then there is no intersection
    if (λ < zero(T)) || (λ > one(T))
      return @IT NoIntersection nothing f
    else
      return @IT CrossingSegmentPlane s(λ) f
    end
  end
end