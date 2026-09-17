# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EdgeRefinement([pred])

Refine polygonal mesh by splitting the edges for which the predicate
`pred` holds true. n-gons with split edges are subdivided into triangles,
and all other n-gons are preserved. Midpoints of split edges are connected
to a new point inside the n-gon when n > 3 or to other midpoints when n == 3.
The default predicate is `e -> true`, meaning all edges are splitted in the
refined mesh.

## Examples

```julia
EdgeRefinement(e -> length(e) > 500u"km")
```
"""
struct EdgeRefinement{F} <: RefinementMethod
  pred::F
end

EdgeRefinement() = EdgeRefinement(e -> true)

function refine(mesh, method::EdgeRefinement)
  assertion(paramdim(mesh) == 2, "EdgeRefinement only defined for surface meshes")

  # convert topology to half-edge structure
  t = convert(HalfEdgeTopology, topology(mesh))

  # original vertices
  vpts = vertices(mesh)

  # midpoints of edges that satisfy the predicate
  mpts = empty(vpts)
  mdict = Dict{Tuple{Int,Int},Int}()
  ∂₁₀ = Boundary{1,0}(t)
  for eind in 1:nfacets(t)
    i, j = ∂₁₀(eind)
    edge = Segment(vpts[i], vpts[j])
    if method.pred(edge)
      push!(mpts, centroid(edge))
      mdict[minmax(i, j)] = length(vpts) + length(mpts)
    end
  end

  # centroids of n-gons with split edges
  cpts = empty(vpts)

  # construct subelements of faces
  ∂₂₀ = Boundary{2,0}(t)
  ngons = Tuple[]
  for eind in 1:nelements(t)
    vinds = ∂₂₀(eind)
    minds = _midpoints(mdict, vinds)
    if all(iszero, minds)
      push!(ngons, vinds)
    elseif length(vinds) == 3
      _subtriangles!(ngons, vinds, minds)
    else
      push!(cpts, coordmean(vpts[i] for i in vinds))
      _subngons!(ngons, vinds, minds, length(vpts) + length(mpts) + length(cpts))
    end
  end

  # new points in refined mesh
  newpoints = [vpts; mpts; cpts]

  # new connectivity in refined mesh
  newconnec = map(connect, ngons)

  SimpleMesh(newpoints, newconnec)
end

# indices of the midpoints inserted on the edges of an element, in the
# order of its vertices, which are zero if the edge is not split. the
# midpoints are shared with the adjacent elements, so the refined mesh
# has no hanging vertices.
function _midpoints(mdict, vinds)
  nv = length(vinds)
  map(1:nv) do i
    u, v = vinds[i], vinds[mod1(i + 1, nv)]
    get(mdict, minmax(u, v), 0)
  end
end

# subdivide the triangle (i, j, k) into two, three or four triangles
# given the midpoints (m1, m2, m3) of its edges (i, j), (j, k) and (k, i)
function _subtriangles!(ngons, (i, j, k), (m1, m2, m3))
  if m1 > 0 && m2 > 0 && m3 > 0
    push!(ngons, (i, m1, m3), (j, m2, m1), (k, m3, m2), (m1, m2, m3))
  elseif m1 > 0 && m2 > 0
    push!(ngons, (k, i, m1), (k, m1, m2), (m2, m1, j))
  elseif m2 > 0 && m3 > 0
    push!(ngons, (i, j, m2), (i, m2, m3), (m3, m2, k))
  elseif m3 > 0 && m1 > 0
    push!(ngons, (j, k, m3), (j, m3, m1), (m1, m3, i))
  elseif m1 > 0
    push!(ngons, (i, m1, k), (m1, j, k))
  elseif m2 > 0
    push!(ngons, (j, m2, i), (m2, k, i))
  else
    push!(ngons, (k, m3, j), (m3, i, j))
  end
end

# subdivide the n-gon into triangles that connect its centroid `c` to
# the vertices and midpoints along its boundary
function _subngons!(ngons, vinds, minds, c)
  ring = Int[]
  for (i, v) in enumerate(vinds)
    push!(ring, v)
    minds[i] > 0 && push!(ring, minds[i])
  end
  nr = length(ring)
  for i in 1:nr
    push!(ngons, (c, ring[i], ring[mod1(i + 1, nr)]))
  end
end
