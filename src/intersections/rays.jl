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
function intersection(f, r1::Ray{N,T}, r2::Ray{N,T}) where {N,T}
  a, b = r1(0), r1(1)
  c, d = r2(0), r2(1)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b, c, d)
  
  # not in same plane or parallel
  if r ≠ rₐ 
    return @IT NoIntersection nothing f #CASE 6
  # collinear
  elseif r == rₐ == 1 
    if direction(r1) ⋅ direction(r2) ≥ 0 # rays aligned in same direction
      if (origin(r1) - origin(r2)) ⋅ direction(r1) ≥ 0 # origin of r1 ∈ r2
        return @IT OverlappingAgreeingRays r1 f # CASE 4: r1
      else
        return @IT OverlappingAgreeingRays r2 f # CASE 4: r2
      end
    else # colliding rays
      if origin(r1) ∉ r2
        return @IT NoIntersection nothing f # CASE 6
      elseif origin(r1) == origin(r2)
        return @IT CornerTouchingRays a f # CASE 3
      else
        return @IT OverlappingOpposingRays Segment(origin(r1), origin(r2)) f # CASE 5
      end
    end
  # in same plane, not parallel
  else
    λ₁ = isapprox(λ₁, zero(T), atol=atol(T)) ? zero(T) : λ₁
    λ₂ = isapprox(λ₂, zero(T), atol=atol(T)) ? zero(T) : λ₂
    if λ₁ < 0 || λ₂ < 0
      return @IT NoIntersection nothing f # CASE 6
    elseif λ₁ == 0
      if λ₂ == 0
        return @IT CornerTouchingRays a f # CASE 3
      else
        return @IT MidTouchingRays a f # CASE 2: origin of r1
      end
    else
      if λ₂ == 0
        return @IT MidTouchingRays c f # CASE 2: origin of r2
      else
        return @IT CrossingRays r1(λ₁) f # CASE 1: equal to r2(λ₂)
      end
    end
  end
end

# compute the intersection of two rays assuming that it is a point
intersectpoint(r1::Ray, r2::Ray) = intersectpoint(Line(r1(0), r1(1)), Line(r2(0), r2(1)))
