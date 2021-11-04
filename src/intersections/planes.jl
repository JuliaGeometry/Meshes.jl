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
    sᵥ = vertices(s)
    p₀ = vertices(p)[1]
    
    n = m.v × m.w
    ln = (sᵥ[2] - sᵥ[1]) ⋅ n
    pₒn = p₀ ⋅ n
  
    if isapprox(ln, zero(T))
        if isapprox(pₒn, ln)
            return ContainedSegmentPlane()
        else
            return NoIntersection()
        end
    else
        λ = ((sᵥ[1] - p₀) ⋅ n) / ln

        λ = isapprox(λ, zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ, one(T), atol=atol(T)) ? one(T) : λ)

        if (λ < zero(T)) || (λ > one(T))
            return IntersectingSegmentPlane(s(λ), λ)
        else
            return NoIntersection()
        end
    end
  end