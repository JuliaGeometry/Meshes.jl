# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of five types:

1. intersect at one inner point (CrossingRaySegment -> Point)
2. intersect at one corner point of segment xor origin of ray (MidTouchingRaySegment -> Point)
3. intersects at one corner point of segment and origin of ray (CornerTouchingRaySegment -> Point)
4. overlap at more than one point (OverlappingRaySegment -> Segment)
5. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, ray::Ray{N,T}, seg::Segment{N,T}) where {N,T}
  a, b = ray(0), ray(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  len₁, len₂ = norm(b - a), length(seg)
  b₀ = a + 1/len₁ * (b - a)
  d₀ = c + 1/len₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NoIntersection nothing f # CASE 5
  # collinear
  elseif r == rₐ == 1
    rc  = sum((c - a) ./ direction(ray))/N
    rd = sum((d - a) ./ direction(ray))/N
    rc = mayberound(rc, zero(T))
    rd = mayberound(rd, zero(T))
    if rc > 0 # c ∈ ray
      if rd ≥ 0
        return @IT OverlappingRaySegment seg f # CASE 4
      else
        return @IT OverlappingRaySegment Segment(origin(ray), c) f # CASE 4
      end
    elseif rc == 0
      if rd > 0
        return @IT OverlappingRaySegment seg f # CASE 4
      else
        return @IT CornerTouchingRaySegment a f # CASE 3
      end
    else # rc < 0
      if rd > 0
        return @IT OverlappingRaySegment (Segment(origin(ray), d)) f # CASE 4
      elseif rd == 0
        return @IT CornerTouchingRaySegment a f # CASE 3
      else
        return @IT NoIntersection nothing f
      end
    end
  # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(mayberound(λ₂, zero(T)), len₂)
    if λ₁ < 0 || (λ₂ < 0 || λ₂ > len₂)
      return @IT NoIntersection nothing f
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == len₂
        return @IT CornerTouchingRaySegment a f # CASE 3
      else
        return @IT MidTouchingRaySegment a f # CASE 2
      end
    else
      if λ₂ == 0 || λ₂ == len₂
        return @IT MidTouchingRaySegment (λ₂ < len₂/2 ? c : d)  f # CASE 2
      else
        return @IT CrossingRaySegment ray(λ₁/len₁) f # CASE 1, equal to seg(λ₂/len₂)
      end
    end
  end
end

