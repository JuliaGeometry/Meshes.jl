# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CatmullClark()

Catmull-Clark refinement of polygonal meshes.

Strictly speaking, the Catmull-Clark algorithm
is used for subdivision surface modeling, not
just mesh refinement. At each step of refinment,
the vertices are adjusted to approximate a smooth
surface.

## References

* Catmull & Clark. 1978. [Recursively generated
  B-spline surfaces on arbitrary topological meshes]
  (https://www.sciencedirect.com/science/article/abs/pii/0010448578901100)
"""
struct CatmullClark <: RefinementMethod end

function refine(mesh, ::CatmullClark)
  # retrieve geometry and topology
  points = vertices(mesh)
  connec = topology(mesh)

  # convert to half-edge structure
  s = convert(HalfEdgeStructure, connec)

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(s)
  epts = map(1:nelements(s)) do elem
    ps = view(points, ∂₂₀(elem))
    cₒ = sum(coordinates, ps) / length(ps)
    Point(cₒ)
  end

  # add midpoints of edges
  ∂₁₂ = Coboundary{1,2}(s)
  ∂₁₀ = Boundary{1,0}(s)
  fpts = map(1:nfacets(s)) do edge
    ps = view(epts, ∂₁₂(edge))
    qs = view(points, ∂₁₀(edge))
    ∑p = sum(coordinates, ps)
    ∑q = sum(coordinates, qs)
    M  = length(ps) + length(qs)
    Point((∑p + ∑q) / M)
  end

  # move original vertices
  ∂₀₂ = Coboundary{0,2}(s)
  ∂₀₀ = Adjacency{0}(s)
  vpts = map(1:nvertices(s)) do u
    # original point
    P = coordinates(points[u])

    # average of centroids
    ps = view(epts, ∂₀₂(u))
    F  = sum(coordinates, ps) / length(ps)

    # average of midpoints
    vs = ∂₀₀(u)
    n  = length(vs)
    R  = sum(vs) do v
      uv = view(points, [u,v])
      sum(coordinates, uv) / 2
    end / n

    Point((F + 2R + (n - 3)P) / n)
  end

  # new points in refined mesh
  newpoints = [vpts; epts; fpts]

  offset₁ = length(vpts)
  offset₂ = offset₁ + length(epts)

  # connect vertices into new quadrangles
  ∂₂₁ = Boundary{2,1}(s)
  newconnec = Connectivity{Quadrangle,4}[]
  for elem in 1:nelements(s)
    verts = CircularVector(∂₂₀(elem))
    edges = CircularVector(∂₂₁(elem))
    for i in 1:length(edges)
      u = elem       + offset₁
      v = edges[i]   + offset₂
      w = verts[i+1]
      z = edges[i+1] + offset₂
      quad = connect((u, v, w, z))
      push!(newconnec, quad)
    end
  end

  SimpleMesh(newpoints, newconnec)
end
