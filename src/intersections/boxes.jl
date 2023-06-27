# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of four types:
# 
# 1. overlap with non-zero measure (Overlapping -> Box)
# 2. intersect at one of the facets (Touching -> Box)
# 3. intersect at corner point (CornerTouching -> Point)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, box₁::Box{Dim,T}, box₂::Box{Dim,T}) where {Dim,T}
  m1, M1 = coordinates.(extrema(box₁))
  m2, M2 = coordinates.(extrema(box₂))

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  # branch on possible configurations
  if u ≺ v
    return @IT Overlapping Box(u, v) f
  elseif u ≻ v
    return @IT NotIntersecting nothing f
  elseif isapprox(u, v, atol=atol(T))
    return @IT CornerTouching u f
  else
    return @IT Touching Box(u, v) f
  end
end
