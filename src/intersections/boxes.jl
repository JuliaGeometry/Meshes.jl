# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of four types:
# 
# 1. overlap with non-zero measure
# 2. intersect at one of the boundaries
# 3. intersect at corner point
# 4. do not overlap nor intersect
function intersection(f, b1::Box{Dim,T}, b2::Box{Dim,T}) where {Dim,T}
  m1, M1 = coordinates.(extrema(b1))
  m2, M2 = coordinates.(extrema(b2))

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  # branch on possible configurations
  if u ≺ v
    return @IT OverlappingBoxes Box(u, v) f
  elseif u ≻ v
    return @IT NoIntersection nothing f
  elseif isapprox(u, v, atol=atol(T))
    return @IT CornerTouchingBoxes u f
  else
    return @IT FaceTouchingBoxes Box(u, v) f
  end
end
