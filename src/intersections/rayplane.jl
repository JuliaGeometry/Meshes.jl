# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
(https://en.wikipedia.org/wiki/Line-plane_intersection)
=#
function Meshes.intersection(f, r::Ray{3,T}, p::Plane{T}) where {T}
    pₒ = coordinates(origin(p))
    n  = normal(p)
    
    # calculate components
    ln = r.v ⋅ n
  
    # if ln is zero, the segment is parallel to the plane
    if isapprox(ln, zero(T), atol=atol(T))
        # if the numerator is zero, the segment is coincident
        if isapprox(coordinates(pₒ - r.p) ⋅ n, zero(T), atol=atol(T))
            return @IT OverlappingRayPlane r f
        else
            return @IT NoIntersection nothing f
        end
    else
        # calculate the ray parameter
        λ = -(n ⋅ coordinates(r.p - pₒ)) / ln

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
end