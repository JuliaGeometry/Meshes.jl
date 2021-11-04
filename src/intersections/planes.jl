# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(s, p)

Compute the intersection of a segment `s` and a plane `p`

## References

* https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
"""
function intersecttype(s::Segment{3,T}, p::Plane{3,T}) where {T}
    # Convert the positions of the segment vertices and the plane origin to coordinates
    sᵥ = coordinates.(vertices(s))
    pₒ = coordinates(origin(p))

    # Get the normal of the plane
    n = normal(p)
    
    # Calculate components
    ln = (sᵥ[2] - sᵥ[1]) ⋅ n
    pₒn = pₒ ⋅ n
  

    if isapprox(ln, zero(T))
        if isapprox(pₒn, ln)
            return ContainedSegmentPlane()
        else
            return NoIntersection()
        end
    else
        # Calculate the segment parameter
        λ = ((sᵥ[1] - pₒ) ⋅ n) / ln

        # If λ is approximately 0 or 1, set as so to prevent any domain errors
        λ = isapprox(λ, zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ, one(T), atol=atol(T)) ? one(T) : λ)

        # If λ is out of bounds for the segment, then there is no intersection
        if (λ < zero(T)) || (λ > one(T))
            return NoIntersection()
        else
            return IntersectingSegmentPlane(s(λ), λ)
        end
    end
  end