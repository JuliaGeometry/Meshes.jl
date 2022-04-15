# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersection(geometry, polygon)

Compute the intersection of a `geometry` and a `polygon`.
"""
function intersection(g, p::Polygon)
  for t in triangulate(p)
    res = intersection(g, t)
    if res.type != NoIntersection
      return res
    end
  end
  return @IT NoIntersection nothing
end