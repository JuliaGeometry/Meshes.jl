# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of six types:

1. intersect at one inner point (CrossingRays -> Point)
2. intersect at origin of one ray (MidTouchingRays -> Point)
3. intersects at origin of both rays (CornerTouchingRays -> Point)
4. overlap with aligned vectors (OverlappingAlignedRays -> Ray)
5. overlap with colliding vectors (OverlappingCollidingRays -> Segment)
6. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, r1::Ray{3,T}, r2::Ray{3,T}) where {T}
  a, b = r1(0), r1(1)
  c, d = r2(0), r2(1)

  if measure(Tetrahedron(a, b, c, d)) > 0 # not in same plane
    return @IT NoIntersection nothing f
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) # parallel
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
      if sum(r1.v .* r2.v) > 0 #aligned rays
        if sum((r1.p - r2.p) .* r1.v) > 0 # origin of r1 ∈ r2
          return @IT OverlappingAlignedRays r1 f
        else
          return @IT OverlappingAlignedRays r2 f
        end
      else #colliding rays
        if r1.p ∉ r2 
          return @IT NoIntersection nothing f
        elseif r1.p == r2.p
          return @IT CornerTouchingRays r1.p f
        else
          return @IT OverlappingCollidingRays Segment(r1.p, r2.p) f
        end
      end
    else # parallel lines, d > 0
      return @IT NoIntersection nothing f
    end
  else # in same plane, not parallel
    v1  = a - b
    v2  = c - d
    v12 = a - c 
    #TODO: find better way to calculate parameters in 3D including third component
    p1 = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    p2 = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    if r1 < 0 || r2 < 0
      return @IT NoIntersection nothing f
    elseif isapprox(p1, zero(T), atol=atol(T))
      if isapprox(p2, zero(T), atol=atol(T)) 
        return @IT CornerTouchingRays a f #CornerTouchingRays # a == c
      else
        return @IT MidTouchingRays a f #MidTouchingRays
      end
    else
      if isapprox(p2, zero(T), atol=atol(T))
        return @IT MidTouchingRays c f# MidTouchingRays #origin of r2
      else
        return @IT CrossingRaySegment (a-p1*v1) f #CrossingRays # returns point # equal to c - s * v2
      end
    end
  end
end
#=
The intersection type can be one of six types:

1. intersect at one inner point (CrossingRays -> Point)
2. intersect at origin of one ray (MidTouchingRays -> Point)
3. intersects at origin of both rays (CornerTouchingRays -> Point)
4. overlap with aligned vectors (OverlappingAlignedRays -> Ray)
5. overlap with colliding vectors (OverlappingCollidingRays -> Segment)
6. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, r1::Ray{2,T}, r2::Ray{2,T}) where {T}
  a, b = r1(0), r1(1)
  c, d = r2(0), r2(1)
  
  if isapprox(abs((b - a) × (c - d)), zero(T), atol=atol(T)^2) # parallel
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
      if sum(r1.v .* r2.v) > 0 #aligned rays
        if sum((r1.p - r2.p) .* r1.v) > 0 # origin of r1 ∈ r2
          return @IT OverlappingAlignedRays r1 f
        else 
          return @IT OverlappingAlignedRays r2 f
        end
      else #colliding rays
        if r1.p ∉ r2 
          return @IT NoIntersection nothing f 
        elseif r1.p == r2.p
          return @IT CornerTouchingRays r1.p f
        else
          return @IT OverlappingCollidingRays Segment(r1.p, r2.p) f
        end
      end
    else # parallel lines, d > 0
      return @IT NoIntersection nothing f
    end
  else
    v1  = a - b
    v2  = c - d
    v12 = a - c 
    #TODO: find better way to calculate parameters in 3D including third component
    p1 = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    p2 = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    if r1 < 0 || r2 < 0
      return @IT NoIntersection nothing f
    elseif isapprox(p1, zero(T), atol=atol(T))
      if isapprox(p2, zero(T), atol=atol(T)) 
        return @IT CornerTouchingRays a f #CornerTouchingRays # a == c
      else
        return @IT MidTouchingRays a f #MidTouchingRays
      end
    else
      if isapprox(p2, zero(T), atol=atol(T))
        return @IT MidTouchingRays c f# MidTouchingRays #origin of r2
      else
        return @IT CrossingRaySegment (a-p1*v1) f #CrossingRays # returns point # equal to c - s * v2
      end
    end
  end
end

# compute the intersection of two rays assuming that it is a point
function intersectpoint(r1::Ray, r2::Ray)
  intersectpoint(Line(r1(0), r1(1)), Line(r2(0), r2(1)))
end
