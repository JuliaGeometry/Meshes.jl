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
function intersection(f, r1::Ray{N,T}, s1::Segment{N,T}) where {N,T}
  a, b = r1(0), r1(1)
  c, d = s1(0), s1(1)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b, c, d)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NoIntersection nothing f # CASE 5
  # collinear
  elseif r == rₐ == 1
    rc  = sum((c - a) ./ direction(r1))/N
    rd = sum((d - a) ./ direction(r1))/N
    rc = isapprox(rc, zero(T), atol=atol(T)) ? zero(T) : rc
    rd = isapprox(rd, zero(T), atol=atol(T)) ? zero(T) : rd
    if rc > 0 # c ∈ r1
      if rd ≥ 0
        return @IT OverlappingRaySegment s1 f # CASE 4
      else
        return @IT OverlappingRaySegment Segment(origin(r1), c) f # CASE 4
      end
    elseif rc == 0
      if rd > 0
        return @IT OverlappingRaySegment s1 f # CASE 4
      else
        return @IT CornerTouchingRaySegment a f # CASE 3
      end
    else # rc < 0
      if rd > 0
        return @IT OverlappingRaySegment (Segment(origin(r1), d)) f # CASE 4
      elseif rd == 0
        return @IT CornerTouchingRaySegment a f # CASE 3
      else
        return @IT NoIntersection nothing f
      end
    end
  # in same plane, not parallel
  else
    λ₁ = isapprox(λ₁, zero(T), atol=atol(T)) ? zero(T) : λ₁
    λ₂ = isapprox(λ₂, zero(T), atol=atol(T)) ? zero(T) : λ₂
    λ₂ = isapprox(λ₂, one(T), atol=atol(T)) ? one(T) : λ₂
    if λ₁ < 0 || (λ₂ < 0 || λ₂ > 1)
      return @IT NoIntersection nothing f
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == 1
        return @IT CornerTouchingRaySegment a f # CASE 3
      else
        return @IT MidTouchingRaySegment a f # CASE 2
      end
    else
      if λ₂ == 0 || λ₂ == 1
        return @IT MidTouchingRaySegment (λ₂ < 0.5 ? c : d)  f # CASE 2
      else
        return @IT CrossingRaySegment r1(λ₁) f # CASE 1, equal to s1(λ₂)
      end
    end
  end
end

# for 2D and 3D use lines.jl implementation
# NOTE: no check whether resulting point is in ray and segment
intersectpoint(r1::Ray, s1::Segment) = intersectpoint(Line(r1(0), r1(1)), Line(s1(0), s1(1)))
