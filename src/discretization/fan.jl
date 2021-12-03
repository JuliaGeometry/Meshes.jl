# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FanTriangulation()

The fan triangulation algorithm fox convex polygons.
See [https://en.wikipedia.org/wiki/Fan_triangulation]
(https://en.wikipedia.org/wiki/Fan_triangulation).
"""
struct FanTriangulation <: DiscretizationMethod end

function discretize(𝒫::Chain, ::FanTriangulation)
  points = vertices(𝒫)
  connec = [connect((1,i,i+1)) for i in 2:nvertices(𝒫)-1]
  SimpleMesh(collect(points), connec)
end