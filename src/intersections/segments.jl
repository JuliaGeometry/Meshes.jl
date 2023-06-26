# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of five types:
#
# 1. intersect at one inner point (CrossingSegments -> Point)
# 2. intersect at one endpoint of one segment (MidTouchingSegments -> Point)
# 3. intersect at one endpoint of both segments (CornerTouchingSegments -> Point)
# 4. overlap of segments (OverlappingSegments -> Segments)
# 5. do not overlap nor intersect (NoIntersection)
function intersection(f, seg₁::Segment{N,T}, seg₂::Segment{N,T}) where {N,T}
  a, b = seg₁(0), seg₁(1)
  c, d = seg₂(0), seg₂(1)

  l₁, l₂ = length(seg₁), length(seg₂)
  b₀ = a + 1 / l₁ * (b - a) # corresponds to seg₁(1/length)
  d₀ = c + 1 / l₂ * (d - c)

  # arc length parameters λ₁ ∈ [0, l₁], λ₂ ∈ [0, l₂]: 
  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f #CASE 5
  elseif r == rₐ == 1 # collinear
    # find parameters λc and λd for points c and d in seg₁
    # use dimension with largest vector component to avoid division by zero
    v = b₀ - a
    i = argmax(abs.(v))
    λc, λd = ((c - a)[i], (d - a)[i]) ./ v[i]
    λc = mayberound(mayberound(λc, zero(T)), l₁)
    λd = mayberound(mayberound(λd, zero(T)), l₁)
    if (λc > l₁ && λd > l₁) || (λc < 0 && λd < 0)
      return @IT NoIntersection nothing f # CASE 5
    elseif (λc == 0 && λd < 0) || (λd == 0 && λc < 0)
      return @IT CornerTouchingSegments a f # CASE 3
    elseif (λc == l₁ && λd > l₁) || (λd == l₁ && λc > l₁)
      return @IT CornerTouchingSegments b f # CASE 3
    else
      params = sort([0, 1, λc / l₁, λd / l₁])
      p₁ = seg₁(params[2])
      p₂ = seg₁(params[3])
      return @IT OverlappingSegments Segment(p₁, p₂) f # CASE 4
    end
  else # in same plane, not parallel
    λ₁ = mayberound(mayberound(λ₁, zero(T)), l₁)
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₁ < 0 || λ₂ < 0 || λ₁ > l₁ || λ₂ > l₂
      return @IT NoIntersection nothing f # CASE 5
    # 8 cases remain
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouchingSegments a f # CASE 3
      else
        return @IT MidTouchingSegments a f # CASE 2
      end
    elseif λ₁ == l₁
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouchingSegments b f # CASE 3
      else
        return @IT MidTouchingSegments b f # CASE 2
      end
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT MidTouchingSegments (λ₂ == 0 ? c : d) f # CASE 2
    else
      return @IT CrossingSegments seg₁(λ₁ / l₁) f # CASE 1: equal to seg₂(λ₂/l₂)
    end
  end
end
