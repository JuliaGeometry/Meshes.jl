# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of six types:

1. intersect at one inner point (CrossingRays -> Point)
2. intersect at origin of one ray (MidTouchingRays -> Point)
3. intersect at origin of both rays (CornerTouchingRays -> Point)
4. overlap with aligned vectors (OverlappingAgreeingRays -> Ray)
5. overlap with colliding vectors (OverlappingOpposingRays -> Segment)
6. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, r1::Ray{3,T}, r2::Ray{3,T}) where {T}
  a, b = r1(0), r1(1)
  c, d = r2(0), r2(1)

  if measure(Tetrahedron(a, b, c, d)) > 0 # not in same plane
    return @IT NoIntersection nothing f #CASE 6
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) # parallel
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
      if r1.v ⋅ r2.v ≥ 0 #rays aligned in same direction
        if (r1.p - r2.p) ⋅ r1.v ≥ 0 # origin of r1 ∈ r2
          return @IT OverlappingAgreeingRays r1 f #CASE 4: r1
        else
          return @IT OverlappingAgreeingRays r2 f #CASE 4: r2
        end
      else #colliding rays
        if r1.p ∉ r2
          return @IT NoIntersection nothing f #CASE 6
        elseif r1.p == r2.p
          return @IT CornerTouchingRays a f #CASE 3
        else
          return @IT OverlappingOpposingRays Segment(r1.p, r2.p) f #CASE 5
        end
      end
    else # parallel lines, d > 0
      return @IT NoIntersection nothing f #CASE 6
    end
  else # in same plane, not parallel
    p1, p2 = intersectparameters(a,b,c,d)
    if p1 < 0 || p2 < 0
      return @IT NoIntersection nothing f #CASE 6
    elseif isapprox(p1, zero(T), atol=atol(T))
      if isapprox(p2, zero(T), atol=atol(T))
        return @IT CornerTouchingRays a f #CASE 3
      else
        return @IT MidTouchingRays a f #CASE 2: origin of r1
      end
    else
      if isapprox(p2, zero(T), atol=atol(T))
        return @IT MidTouchingRays c f # CASE 2: origin of r2
      else
        return @IT CrossingRays (a-p1*v1) f #CASE 1: equal to c - p2 * v2
      end
    end
  end
end

#=
The intersection type can be one of six types:

1. intersect at one inner point (CrossingRays -> Point)
2. intersect at origin of one ray (MidTouchingRays -> Point)
3. intersect at origin of both rays (CornerTouchingRays -> Point)
4. overlap with aligned vectors (OverlappingAgreeingRays -> Ray)
5. overlap with colliding vectors (OverlappingOpposingRays -> Segment)
6. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, r1::Ray{2,T}, r2::Ray{2,T}) where {T}
  a, b = r1(0), r1(1)
  c, d = r2(0), r2(1)

  if isapprox(abs((b - a) × (c - d)), zero(T), atol=atol(T)^2) # parallel
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
      if r1.v ⋅ r2.v ≥ 0 #aligned rays
        if (r1.p - r2.p) ⋅ r1.v ≥ 0 # origin of r1 ∈ r2
          return @IT OverlappingAgreeingRays r1 f #CASE 4: r1
        else
          return @IT OverlappingAgreeingRays r2 f #CASE 4: r2
        end
      else #colliding rays
        if r1.p ∉ r2 
          return @IT NoIntersection nothing f #CASE 6
        elseif r1.p == r2.p
          return @IT CornerTouchingRays a f #CASE 3
        else
          return @IT OverlappingOpposingRays Segment(r1.p, r2.p) f #CASE 5
        end
      end
    else # parallel lines, d > 0
      return @IT NoIntersection nothing f #CASE 6
    end
  else
    p1, p2 = intersectparameters(a,b,c,d)
    if p1 < 0 || p2 < 0
      return @IT NoIntersection nothing f
    elseif isapprox(p1, zero(T), atol=atol(T))
      if isapprox(p2, zero(T), atol=atol(T))
        return @IT CornerTouchingRays a f #CASE 3
      else
        return @IT MidTouchingRays a f #CASE 2: origin of r1
      end
    else
      if isapprox(p2, zero(T), atol=atol(T))
        return @IT MidTouchingRays c f #CASE 2: origin of r2
      else
        return @IT CrossingRays (a-p1*v1) f #CASE 1: equal to c - p2 * v2
      end
    end
  end
end

# compute the intersection of two rays assuming that it is a point
intersectpoint(r1::Ray, r2::Ray) = intersectpoint(Line(r1.p, r1.p + r1.v), Line(r2.p, r2.p + r2.v))
