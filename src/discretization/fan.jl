# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FanTriangulation()

The fan triangulation algorithm for convex polygons.
See [https://en.wikipedia.org/wiki/Fan_triangulation]
(https://en.wikipedia.org/wiki/Fan_triangulation).
"""
struct FanTriangulation <: BoundaryDiscretizationMethod end

discretizewithin(chain::Chain{2}, ::FanTriangulation) = fan(chain)

discretizewithin(chain::Chain{3}, ::FanTriangulation) = fan(chain)

function fan(chain::Chain)
  points = vertices(chain)
  connec = [connect((1, i, i + 1)) for i in 2:nvertices(chain)-1]
  SimpleMesh(collect(points), connec)
end
