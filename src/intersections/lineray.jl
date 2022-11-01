# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of four types:

1. intersect at one inner point (CrossingLineRay -> Point)
2. intersect at origin of ray (TouchingLineRay -> Point)
3. overlap of line and ray (OverlappingLineRay -> Ray)
4. do not overlap nor intersect (NoIntersection)
=#
function intersection(f, l1::Line{N,T}, r1::Ray{N,T}) where {N,T}
  a, b = l1(0), l1(1)
  c, d = r1(0), r1(1)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b, c, d)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT OverlappingLineRay r1 f # CASE 3
  else # in same plane, not parallel
    λ₂ = isapprox(λ₂, 0, atol=atol(T)) ? 0 : λ₂
    if λ₂ > 0
      return @IT CrossingLineRay l1(λ₁) f # CASE 1
    elseif λ₂ == 0
      return @IT TouchingLineRay origin(r1) f # CASE 2
    else
      return @IT NoIntersection nothing f # CASE 4
    end
  end
end
