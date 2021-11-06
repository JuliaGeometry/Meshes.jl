# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(segment, plane)

Compute the intersection of a segment `s` and a plane `p`
See https://en.wikipedia.org/wiki/Line-plane_intersection
"""
function intersecttype(s::Segment{3,T}, p::Plane{3,T}) where {T}
    # convert the positions of the segment vertices and the plane origin to coordinates
    sᵥ = coordinates.(vertices(s))
    pₒ = coordinates(origin(p))

    # get the normal of the plane
    n = normal(p)
    
    # calculate components
    ln = (sᵥ[2] - sᵥ[1]) ⋅ n
    pₒn = pₒ ⋅ n
  
    # if ln is zero, the segment is parallel to the plane
    if isapprox(ln, zero(T), atol=atol(T))
        # if the numerator is zero, the segment is coincident
        if isapprox((pₒ - sᵥ[1]) ⋅ n, zero(T), atol=atol(T))
            return OverlappingSegmentPlane(s)
        else
            return NoIntersection()
        end
    else
        # calculate the segment parameter
        λ = ((pₒ - sᵥ[1]) ⋅ n) / ln

        # if λ is approximately 0 or 1, set as so to prevent any domain errors
        if isapprox(λ, zero(T), atol=atol(T))
            return TouchingSegmentPlane(s(zero(T)))
        elseif isapprox(λ, one(T), atol=atol(T))
            return TouchingSegmentPlane(s(one(T)))
        end

        # if λ is out of bounds for the segment, then there is no intersection
        if (λ < zero(T)) || (λ > one(T))
            return NoIntersection()
        else
            return CrossingSegmentPlane(s(λ))
        end
    end
  end