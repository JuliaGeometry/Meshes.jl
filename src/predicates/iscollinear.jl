# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    iscollinear(A, B, C)

Tells whether or not the points `A`, `B` and `C` are collinear.
"""
function iscollinear(A::Point, B::Point, C::Point)
  # points A, B, C are collinear if and only if the
  # cross-products for segments AB and AC with respect
  # to all possible pairs of coordinates are zero
  Dim = embeddim(A)
  AB, AC = B - A, C - A
  result = true
  for i in 1:Dim, j in (i + 1):Dim
    u = Vec(AB[i], AB[j])
    v = Vec(AC[i], AC[j])
    if !isapproxzero(u × v)
      result = false
      break
    end
  end
  result
end
