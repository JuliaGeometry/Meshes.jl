# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EdgeRefinement(pred)

Refinement of a mesh by splitting the edges for which the predicate
`pred` holds true. The mesh is triangulated if needed, triangles with
split edges are subdivided into two, three or four triangles, and all
other triangles are preserved.

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

  # triangulate mesh if necessary
  tmesh = eltype(mesh) <: Triangle ? mesh : simplexify(mesh)

  # convert topology to half-edge structure
  t = convert(HalfEdgeTopology, topology(tmesh))

  # original vertices
  vpts = vertices(tmesh)

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

  # new points in refined mesh
  newpoints = [vpts; mpts]

  # construct subtriangles of faces
  ∂₂₀ = Boundary{2,0}(t)
  triangles = Tuple{Int,Int,Int}[]
  for tind in 1:nelements(t)
    i, j, k = ∂₂₀(tind)
    m1 = get(mdict, minmax(i, j), 0)
    m2 = get(mdict, minmax(j, k), 0)
    m3 = get(mdict, minmax(k, i), 0)
    _subtriangles!(triangles, (i, j, k), (m1, m2, m3))
  end

  SimpleMesh(newpoints, map(connect, triangles))
end

# subdivide the triangle (i, j, k) given the midpoints (m1, m2, m3)
# of its edges (i, j), (j, k) and (k, i), which are zero if the edge
# is not split. the midpoints are shared with the adjacent triangles,
# so the resulting mesh has no hanging vertices.
function _subtriangles!(triangles, (i, j, k), (m1, m2, m3))
  if m1 > 0 && m2 > 0 && m3 > 0
    push!(triangles, (i, m1, m3), (j, m2, m1), (k, m3, m2), (m1, m2, m3))
  elseif m1 > 0 && m2 > 0
    push!(triangles, (k, i, m1), (k, m1, m2), (m2, m1, j))
  elseif m2 > 0 && m3 > 0
    push!(triangles, (i, j, m2), (i, m2, m3), (m3, m2, k))
  elseif m3 > 0 && m1 > 0
    push!(triangles, (j, k, m3), (j, m3, m1), (m1, m3, i))
  elseif m1 > 0
    push!(triangles, (i, m1, k), (m1, j, k))
  elseif m2 > 0
    push!(triangles, (j, m2, i), (m2, k, i))
  elseif m3 > 0
    push!(triangles, (k, m3, j), (m3, i, j))
  else
    push!(triangles, (i, j, k))
  end
end
