# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ManualSimplexification()

Simplexify geometries manually using indices of vertices.
"""
struct ManualSimplexification <: DiscretizationMethod end

discretize(box::Box{𝔼{1}}, ::ManualSimplexification) = SimpleMesh(collect(extrema(box)), GridTopology(1))

function discretize(box::Box{𝔼{2}}, ::ManualSimplexification)
  indices = [(1, 2, 3), (1, 3, 4)]
  SimpleMesh(pointify(box), connect.(indices, Triangle))
end

function discretize(box::Box{𝔼{3}}, ::ManualSimplexification)
  indices = [(1, 5, 6, 8), (1, 3, 4, 8), (1, 3, 6, 8), (1, 2, 3, 6), (3, 6, 7, 8)]
  SimpleMesh(pointify(box), connect.(indices, Tetrahedron))
end

function discretize(box::Box{🌐}, ::ManualSimplexification)
  indices = [(1, 2, 3), (1, 3, 4)]
  SimpleMesh(pointify(box), connect.(indices, Triangle))
end

function discretize(hexa::Hexahedron, ::ManualSimplexification)
  indices = [(1, 5, 6, 8), (1, 3, 4, 8), (1, 3, 6, 8), (1, 2, 3, 6), (3, 6, 7, 8)]
  SimpleMesh(pointify(hexa), connect.(indices, Tetrahedron))
end

function discretize(pyramid::Pyramid, ::ManualSimplexification)
  indices = [(1, 2, 4, 5), (3, 4, 2, 5)]
  SimpleMesh(pointify(pyramid), connect.(indices, Tetrahedron))
end

function discretize(wedge::Wedge, ::ManualSimplexification)
  indices = [(1, 2, 3, 4), (4, 5, 6, 2), (4, 5, 6, 3)]
  SimpleMesh(pointify(wedge), connect.(indices, Tetrahedron))
end
