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
function intersection(f, ray::Ray{N,T}, line::Line{N,T}) where {N,T}
  a, b = ray(0), ray(1)
  c, d = line(0), line(1)

  # rescaling of point b necessary to gain a parameter λ₁ representing the arc length
  len₁ = norm(b - a)
  b₀ = a + 1/len₁ * (b - a)

  λ₁, _, r, rₐ = intersectparameters(a, b₀, c, d)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT OverlappingRayLine ray f # CASE 3
  else # in same plane, not parallel
    λ₁ = mayberound(λ₁, zero(T))
    if λ₁ > 0
      return @IT CrossingRayLine ray(λ₁/len₁) f # CASE 1
    elseif λ₁ == 0
      return @IT TouchingRayLine origin(ray) f # CASE 2
    else
      return @IT NoIntersection nothing f # CASE 4
    end
  end
end
