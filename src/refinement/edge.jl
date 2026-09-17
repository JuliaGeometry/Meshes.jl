# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EdgeRefinement(pred)

Refinement of a mesh by splitting the edges for which the predicate
`pred` holds true. Elements with split edges are subdivided into
triangles, and all other elements are preserved.

## Examples

```julia
EdgeRefinement(e -> measure(e) > 500u"km")
```
"""
struct EdgeRefinement{F} <: RefinementMethod
  pred::F
end

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
  connec = Connectivity[]
  for eind in 1:nelements(t)
    verts = ∂₂₀(eind)
    mids = _midpoints(mdict, verts)
    if all(iszero, mids)
      push!(connec, element(t, eind))
    elseif length(verts) == 3
      _subtriangles!(connec, verts, mids)
    else
      push!(cpts, coordmean(vpts[i] for i in verts))
      _subngons!(connec, verts, mids, length(vpts) + length(mpts) + length(cpts))
    end
  end

  # new points in refined mesh
  newpoints = [vpts; mpts; cpts]

  SimpleMesh(newpoints, map(identity, connec))
end

# indices of the midpoints inserted on the edges of an element, in the
# order of its vertices, which are zero if the edge is not split. the
# midpoints are shared with the adjacent elements, so the refined mesh
# has no hanging vertices.
function _midpoints(mdict, verts)
  nv = length(verts)
  map(1:nv) do i
    u, v = verts[i], verts[mod1(i + 1, nv)]
    get(mdict, minmax(u, v), 0)
  end
end

# subdivide the triangle (i, j, k) into two, three or four triangles
# given the midpoints (m1, m2, m3) of its edges (i, j), (j, k) and (k, i)
function _subtriangles!(connec, (i, j, k), (m1, m2, m3))
  if m1 > 0 && m2 > 0 && m3 > 0
    push!(connec, connect((i, m1, m3)), connect((j, m2, m1)), connect((k, m3, m2)), connect((m1, m2, m3)))
  elseif m1 > 0 && m2 > 0
    push!(connec, connect((k, i, m1)), connect((k, m1, m2)), connect((m2, m1, j)))
  elseif m2 > 0 && m3 > 0
    push!(connec, connect((i, j, m2)), connect((i, m2, m3)), connect((m3, m2, k)))
  elseif m3 > 0 && m1 > 0
    push!(connec, connect((j, k, m3)), connect((j, m3, m1)), connect((m1, m3, i)))
  elseif m1 > 0
    push!(connec, connect((i, m1, k)), connect((m1, j, k)))
  elseif m2 > 0
    push!(connec, connect((j, m2, i)), connect((m2, k, i)))
  else
    push!(connec, connect((k, m3, j)), connect((m3, i, j)))
  end
end

# subdivide the n-gon into triangles that connect its centroid `c` to
# the vertices and midpoints along its boundary
function _subngons!(connec, verts, mids, c)
  ring = Int[]
  for (i, v) in enumerate(verts)
    push!(ring, v)
    mids[i] > 0 && push!(ring, mids[i])
  end
  nr = length(ring)
  for i in 1:nr
    push!(connec, connect((c, ring[i], ring[mod1(i + 1, nr)])))
  end
end
