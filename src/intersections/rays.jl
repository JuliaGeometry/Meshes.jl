# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersecttype(r::Ray{Dim,T}, tri::Triangle{Dim,T}) where {Dim, T}
    # get the triangle's edges and it's normal
    tri_v = vertices(tri)
    u = tri_v[2] - tri_v[1]
    v = tri_v[3] - tri_v[1]
    tri_n = normal(tri)

    # calculate the numerator and denominator to determine the intersection
    # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    n0 = tri_v[1] - r.p 
    a = n0 ⋅ tri_n
    b = r.v ⋅ tri_n

    # if the denominator is 0, then the ray is parallel to the tri
    if isapprox(b, zero(T))
        return NonIntersectingRayTri()
    end

    # find ray parameter at intersection with triangle
    t = a / b

    # if t is less than zero, the triangle is "behind" the ray
    # and doesn't intersect
    if t < 0.0                    
        return NonIntersectingRayTri()                  
    end 
    
    # evaluate the intersection point
    p = r.p + t * r.v

    # if the point is within the plane AND the triangle then it intersects
    if p ∈ tri
        return IntersectingRayTri(p, t)
    else
        return NonIntersectingRayTri()
    end
end
