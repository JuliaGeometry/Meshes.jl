# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ManualSimplexification()

Simplexify convex geometries manually using indices of vertices.
"""
struct ManualSimplexification <: DiscretizationMethod end

discretize(geom::Geometry, ::ManualSimplexification) = SimpleMesh(pointify(geom), _manualconnec(geom))

function discretize(chain::Chain, ::ManualSimplexification)
  np = nvertices(chain) + isclosed(chain)
  ip = isperiodic(chain)

  points = collect(eachvertex(chain))
  topo = GridTopology((np - 1,), ip)

  SimpleMesh(points, topo)
end

function _manualconnec(::Box{𝔼{1}})
  [connect((1, 2), Segment)]
end

function _manualconnec(::Box{𝔼{2}})
  [connect((1, 2, 3), Triangle), connect((1, 3, 4), Triangle)]
end

function _manualconnec(::Box{𝔼{3}})
  [
    connect((1, 5, 6, 8), Tetrahedron),
    connect((1, 3, 4, 8), Tetrahedron),
    connect((1, 3, 6, 8), Tetrahedron),
    connect((1, 2, 3, 6), Tetrahedron),
    connect((3, 6, 7, 8), Tetrahedron)
  ]
end

function _manualconnec(::Box{🌐})
  [connect((1, 2, 3), Triangle), connect((1, 3, 4), Triangle)]
end

function _manualconnec(::Triangle)
  [connect((1, 2, 3), Triangle)]
end

function _manualconnec(::Hexahedron)
  [
    connect((1, 5, 6, 8), Tetrahedron),
    connect((1, 3, 4, 8), Tetrahedron),
    connect((1, 3, 6, 8), Tetrahedron),
    connect((1, 2, 3, 6), Tetrahedron),
    connect((3, 6, 7, 8), Tetrahedron)
  ]
end

function _manualconnec(::Pyramid)
  [connect((1, 2, 4, 5), Tetrahedron), connect((3, 4, 2, 5), Tetrahedron)]
end

function _manualconnec(::Wedge)
  [connect((1, 2, 3, 4), Tetrahedron), connect((4, 5, 6, 2), Tetrahedron), connect((4, 5, 6, 3), Tetrahedron)]
end
