# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of four types:

1. intersect at one inner point (CrossingRayLine -> Point)
2. intersect at origin of ray (TouchingRayLine -> Point)
3. overlap of line and ray (OverlappingRayLine -> Ray)
4. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, r1::Ray{N,T}, l1::Line{N,T}) where {N,T}
  a, b = r1(0), r1(1)
  c, d = l1(0), l1(1)

  # use a normalized point for ray parameter λ₁
  len1 = norm(b - a)
  b₀ = a + 1/len1 * (b - a)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT OverlappingRayLine r1 f # CASE 3
  else # in same plane, not parallel
    λ₁ = mayberound(λ₁, zero(T))
    if λ₁ > 0
      return @IT CrossingRayLine r1(λ₁/len1) f # CASE 1
    elseif λ₁ == 0
      return @IT TouchingRayLine origin(r1) f # CASE 2
    else
      return @IT NoIntersection nothing f # CASE 4
    end
  end
end
