# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of three types:
# 
# 1. intersect at one point
# 2. overlap at more than one point
# 3. do not overlap nor intersect
function intersection(f, line₁::Line, line₂::Line)
  a, b = line₁(0), line₁(1)
  c, d = line₂(0), line₂(1)

  λ₁, _, r, rₐ = intersectparameters(a, b, c, d)

  if r == rₐ == 2
    return @IT Crossing (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT Overlapping line₁ f
  else
    return @IT NotIntersecting nothing f
  end
end
