# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of six types:
1. intersect at one inner point (CrossingLineSegment -> Point)
2. intersect at an end point of segment (TouchingLineSegment -> Point)
3. overlap of line and segment (OverlappingLineSegment -> Segment)
4. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, line::Line{N,T}, seg::Segment{N,T}) where {N,T}
  a, b = line(0), line(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameter λ₂ corresponding to arc lengths
  l₂ = length(seg)
  d₀ = c + 1 / l₂ * (d - c)

  _, λ₂, r, rₐ = intersectparameters(a, b, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NoIntersection nothing f # CASE 4
  # collinear
  elseif r == rₐ == 1
    return @IT OverlappingLineSegment seg f # CASE 3
  # in same plane, not parallel
  else
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₂ > 0 && λ₂ < l₂
      return @IT CrossingLineSegment seg(λ₂ / l₂) f # CASE 1, equal to line(λ₁)
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT TouchingLineSegment ((λ₂ == 0) ? c : d) f # CASE 2
    else
      return @IT NoIntersection nothing f # CASE 4
    end
  end
end
