# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CatmullClarkRefinement()

Catmull-Clark refinement of polygonal meshes.

Strictly speaking, the Catmull-Clark algorithm
is used for subdivision surface modeling, not
just mesh refinement. At each step of refinement,
the vertices are adjusted to approximate a smooth
surface.

## References

* Catmull & Clark. 1978. [Recursively generated
  B-spline surfaces on arbitrary topological meshes]
  (https://www.sciencedirect.com/science/article/abs/pii/0010448578901100)
"""
struct CatmullClarkRefinement <: RefinementMethod end

function refine(mesh, ::CatmullClarkRefinement)
  # retrieve geometry and topology
  points = vertices(mesh)
  connec = topology(mesh)

  # convert to half-edge structure
  t = convert(HalfEdgeTopology, connec)

  # add centroids of elements
  ∂₂₀ = Boundary{2,0}(t)
  epts = map(1:nelements(t)) do elem
    is = ∂₂₀(elem)
    cₒ = sum(j -> to(points[j]), is) / length(is)
    withcrs(mesh, cₒ)
  end

  # add midpoints of edges
  ∂₁₂ = Coboundary{1,2}(t)
  ∂₁₀ = Boundary{1,0}(t)
  fpts = map(1:nfacets(t)) do edge
    is = ∂₁₂(edge)
    js = ∂₁₀(edge)
    ∑p = sum(i -> to(epts[i]), is)
    ∑q = sum(j -> to(points[j]), js)
    M = length(is) + length(js)
    withcrs(mesh, (∑p + ∑q) / M)
  end

  # move original vertices
  ∂₀₂ = Coboundary{0,2}(t)
  ∂₀₀ = Adjacency{0}(t)
  vpts = map(1:nvertices(t)) do u
    # original point
    P = to(points[u])

    # average of centroids
    is = ∂₀₂(u)
    F = sum(i -> to(epts[i]), is) / length(is)

    # average of midpoints
    vs = ∂₀₀(u)
    n = length(vs)
    R = sum(v -> to(points[u]) + to(points[v]), vs) / 2n

    withcrs(mesh, (F + 2R + (n - 3)P) / n)
  end

  # new points in refined mesh
  newpoints = [vpts; epts; fpts]

  offset₁ = length(vpts)
  offset₂ = offset₁ + length(epts)

  # connect vertices into new quadrangles
  ∂₂₁ = Boundary{2,1}(t)
  newconnec = Connectivity{Quadrangle,4}[]
  for elem in 1:nelements(t)
    verts = ∂₂₀(elem)
    edges = ∂₂₁(elem)
    nv = length(verts)
    ne = length(edges)
    for i in 1:ne
      u = elem + offset₁
      v = edges[mod1(i, ne)] + offset₂
      w = verts[mod1(i + 1, nv)]
      z = edges[mod1(i + 1, ne)] + offset₂
      quad = connect((u, v, w, z))
      push!(newconnec, quad)
    end
  end

  SimpleMesh(newpoints, newconnec)
end
