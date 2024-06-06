# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    QuadRefinement()

Refinement of polygonal meshes into quadrangles.
A n-gon is subdivided into n quadrangles.
"""
struct QuadRefinement <: RefinementMethod end

function refine(mesh, ::QuadRefinement)
  # retrieve geometry and topology
  points = vertices(mesh)
  connec = topology(mesh)

  # convert to half-edge structure
  t = convert(HalfEdgeTopology, connec)

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(t)
  epts = map(1:nelements(t)) do elem
    ps = view(points, ∂₂₀(elem))
    cₒ = sum(to, ps) / length(ps)
    withdatum(mesh, cₒ)
  end

  # add midpoints of edges
  ∂₁₀ = Boundary{1,0}(t)
  fpts = map(1:nfacets(t)) do edge
    ps = view(points, ∂₁₀(edge))
    cₒ = sum(to, ps) / length(ps)
    withdatum(mesh, cₒ)
  end

  # original vertices
  vpts = points

  # new points in refined mesh
  newpoints = [vpts; epts; fpts]

  offset₁ = length(vpts)
  offset₂ = offset₁ + length(epts)

  # connect vertices into new quadrangles
  ∂₂₁ = Boundary{2,1}(t)
  newconnec = Connectivity{Quadrangle,4}[]
  for elem in 1:nelements(t)
    verts = CircularVector(∂₂₀(elem))
    edges = CircularVector(∂₂₁(elem))
    for i in 1:length(edges)
      u = elem + offset₁
      v = edges[i] + offset₂
      w = verts[i + 1]
      z = edges[i + 1] + offset₂
      quad = connect((u, v, w, z))
      push!(newconnec, quad)
    end
  end

  SimpleMesh(newpoints, newconnec)
end
