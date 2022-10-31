# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------
#=
The intersection type can be one of six types:
1. intersect at one inner point (CrossingSegments -> Point)
2. intersect at one endpoint of one segment (MidTouchingSegments -> Point)
3. intersect at one endpoint of both segments (CornerTouchingSegments -> Point)
4. overlap of segments (OverlappingSegments -> Segments)
5. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, s1::Segment{N,T}, s2::Segment{N,T}) where {N,T}
  a, b = s1(0), s1(1)
  c, d = s2(0), s2(1)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b, c, d)
  # helper function to round values close to x (one or zero)
  mayberound(λ, x) = isapprox(λ, x, atol=atol(T)) ? x : λ

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f #CASE 5
  elseif r == rₐ == 1 # collinear
    # find parameters λc and λd for points c and d in s1
    # use dimension with largest vector component to avoid division by zero
    v = b - a
    i = argmax(abs.(v))
    λc, λd = ((c - a)[i], (d - a)[i]) ./ v[i]
    λc = mayberound(mayberound(λc, 0), 1)
    λd = mayberound(mayberound(λd, 0), 1)
    if (λc > 1 && λd > 1) || (λc < 0 && λd < 0)
      return @IT NoIntersection nothing f # CASE 5
    elseif (λc == 0 && λd < 0) ||  (λd == 0 && λc < 0)
      return @IT CornerTouchingSegments a f # CASE 3
    elseif (λc == 1 && λd > 1) || (λd == 1 && λc > 1)
      return @IT CornerTouchingSegments b f # CASE 3
    else
      parameters = sort([0,1,λc,λd])
      return @IT OverlappingSegments Segment(s1(parameters[2]), s1(parameters[3])) f # CASE 4
    end
  else # in same plane, not parallel
    λ₁ = mayberound(mayberound(λ₁, 0), 1)
    λ₂ = mayberound(mayberound(λ₂, 0), 1)
    if λ₁ < 0 || λ₂ < 0 || λ₁ > 1 || λ₂ > 1
      return @IT NoIntersection nothing f # CASE 5
    # 8 cases remain
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == 1
        return @IT CornerTouchingSegments a f # CASE 3
      else
        return @IT MidTouchingSegments a f # CASE 2
      end
    elseif λ₁ == 1
      if λ₂ == 0 || λ₂ == 1
        return @IT CornerTouchingSegments b f # CASE 3
      else
        return @IT MidTouchingSegments b f # CASE 2
      end
    elseif λ₂ == 0 || λ₂ == 1
      return @IT MidTouchingSegments (λ₂ == 0 ? c : d) f # CASE 2
    else
      return @IT CrossingSegments s1(λ₁) f # CASE 1: equal to s2(λ₂)
    end
  end
end