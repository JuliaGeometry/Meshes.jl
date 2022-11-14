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
function intersection(f, ray₁::Ray{N,T}, ray₂::Ray{N,T}) where {N,T}
  a, b = ray₁(0), ray₁(1)
  c, d = ray₂(0), ray₂(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  len₁, len₂ = norm(b - a), norm(d - c)
  b₀ = a + 1/len₁ * (b - a)
  d₀ = c + 1/len₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)
  
  # not in same plane or parallel
  if r ≠ rₐ 
    return @IT NoIntersection nothing f #CASE 6
  # collinear
  elseif r == rₐ == 1 
    if direction(ray₁) ⋅ direction(ray₂) ≥ 0 # rays aligned in same direction
      if (origin(ray₁) - origin(ray₂)) ⋅ direction(ray₁) ≥ 0 # origin of ray₁ ∈ ray₂
        return @IT OverlappingAgreeingRays ray₁ f # CASE 4: ray₁
      else
        return @IT OverlappingAgreeingRays ray₂ f # CASE 4: ray₂
      end
    else # colliding rays
      if origin(ray₁) ∉ ray₂
        return @IT NoIntersection nothing f # CASE 6
      elseif origin(ray₁) == origin(ray₂)
        return @IT CornerTouchingRays a f # CASE 3
      else
        return @IT OverlappingOpposingRays Segment(origin(ray₁), origin(ray₂)) f # CASE 5
      end
    end
  # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(λ₂, zero(T))
    if λ₁ < 0 || λ₂ < 0
      return @IT NoIntersection nothing f # CASE 6
    elseif λ₁ == 0
      if λ₂ == 0
        return @IT CornerTouchingRays a f # CASE 3
      else
        return @IT MidTouchingRays a f # CASE 2: origin of ray₁
      end
    else
      if λ₂ == 0
        return @IT MidTouchingRays c f # CASE 2: origin of ray₂
      else
        return @IT CrossingRays ray₁(λ₁/len₁) f # CASE 1: equal to ray₂(λ₂/len₂)
      end
    end
  end
end
