# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of four types:
# 
# 1. overlap with non-zero measure (Overlapping -> Box)
# 2. intersect at corner point (CornerTouching -> Point)
# 3. intersect at one of the facets (Touching -> Box)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, box₁::Box{Dim}, box₂::Box{Dim}) where {Dim}
  # retrieve corner points
  m1, M1 = coordinates.(extrema(box₁))
  m2, M2 = coordinates.(extrema(box₂))

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  # auxiliary variables
  δ = v - u
  δ̄ = abs.(δ)
  τ = atol(eltype(δ))

  # branch on possible configurations
  if all(>(τ), δ)
    return @IT Overlapping Box(u, v) f
  elseif all(<(τ), δ̄)
    return @IT CornerTouching u f
  elseif any(<(τ), δ̄) && (δ == δ̄ || δ == -δ̄)
    return @IT Touching (u ⪯ v ? Box(u, v) : Box(v, u)) f
  else
    return @IT NotIntersecting nothing f
  end
end
