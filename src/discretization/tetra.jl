# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Tetrahedralization()

A method to discretize geometries into tetrahedra.
"""
struct Tetrahedralization <: DiscretizationMethod end

function discretize(box::Box, ::Tetrahedralization)
  indices = [(1, 5, 6, 8), (1, 3, 4, 8), (1, 3, 6, 8), (1, 2, 3, 6), (3, 6, 7, 8)]
  SimpleMesh(pointify(box), connect.(indices, Tetrahedron))
end

function discretize(hexa::Hexahedron, ::Tetrahedralization)
  indices = [(1, 5, 6, 8), (1, 3, 4, 8), (1, 3, 6, 8), (1, 2, 3, 6), (3, 6, 7, 8)]
  SimpleMesh(pointify(hexa), connect.(indices, Tetrahedron))
end

function discretize(pyramid::Pyramid, ::Tetrahedralization)
  indices = [(1, 2, 4, 5), (3, 4, 2, 5)]
  SimpleMesh(pointify(pyramid), connect.(indices, Tetrahedron))
end

function discretize(wedge::Wedge, ::Tetrahedralization)
  indices = [(1, 2, 3, 4), (4, 5, 6, 2), (4, 5, 6, 3)]
  SimpleMesh(pointify(wedge), connect.(indices, Tetrahedron))
end
