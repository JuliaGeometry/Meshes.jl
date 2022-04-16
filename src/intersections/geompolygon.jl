# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function intersection(f, g::Geometry, p::Polygon)
  for t in triangulate(p)
    I = intersection(g, t)
    if type(I) != NoIntersection
      return I
    end
  end
  return @IT NoIntersection nothing f
end