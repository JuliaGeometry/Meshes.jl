# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TriRefinement()

Refinement of polygonal meshes into triangles.
A n-gon is subdivided into n-2 triangles.
"""
struct TriRefinement <: RefinementMethod end

function refine(mesh, method::TriRefinement)
  if eltype(mesh) <: Triangle
    # go ahead and refine
    _refine(mesh, method)
  else
    # simplexify non-triangle elements
    simplexify(mesh)
  end
end

function _refine(mesh, ::TriRefinement)
  # retrieve geometry and topology
  points = vertices(mesh)
  connec = topology(mesh)

  # convert to half-edge structure
  t = convert(HalfEdgeTopology, connec)

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(t)
  epts = map(1:nelements(t)) do elem
    ps = view(points, ∂₂₀(elem))
    cₒ = sum(coordinates, ps) / length(ps)
    Point(cₒ)
  end

  # original vertices
  vpts = points

  # new points in refined mesh
  newpoints = [vpts; epts]

  offset = length(vpts)

  # connect vertices into new triangles
  newconnec = Connectivity{Triangle,3}[]
  for elem in 1:nelements(t)
    verts = CircularVector(∂₂₀(elem))
    for i in 1:length(verts)
      u = elem + offset
      v = verts[i]
      w = verts[i + 1]
      tri = connect((u, v, w))
      push!(newconnec, tri)
    end
  end

  return SimpleMesh(newpoints, newconnec)
end
