# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    signarea(p₁, p₂, p₃)

Compute signed area of triangle formed
by points `p₁`, `p₂` and `p₃`.
"""
function signarea(p₁::Point{2}, p₂::Point{2}, p₃::Point{2})
  a = coordinates(p₁)
  b = coordinates(p₂)
  c = coordinates(p₃)
  ((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])) / 2
end
