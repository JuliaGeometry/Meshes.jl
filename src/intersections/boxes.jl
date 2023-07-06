# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, box₁::Box, box₂::Box)
  # retrieve corner points
  m1, M1 = coordinates.(extrema(box₁))
  m2, M2 = coordinates.(extrema(box₂))

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  # branch on possible configurations
  if u ≺ v
    return @IT Overlapping Box(u, v) f
  else
    return @IT NotIntersecting nothing f
  end
end
