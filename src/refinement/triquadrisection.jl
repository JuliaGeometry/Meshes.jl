# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TriQuadrisection()

Refinement of a mesh by preliminarly triangulating it if needed and
then subdividing each triangle into four triangles.

## References

* Charles Loop. 1987. [Smooth subdivision surfaces based on 
  triangles](https://charlesloop.com/thesis.pdf). 
  Master's thesis, University of Utah.
"""
struct TriQuadrisection <: RefinementMethod end

function refine(mesh, ::TriQuadrisection)
  if !(eltype(mesh) <: Triangle)
    tmesh = simplexify(mesh)
  else
    tmesh = mesh
  end

  # retrieve vertices
  points = vertices(tmesh)
  npoints = length(points)

  # convert topology to half-edge structure
  t = convert(HalfEdgeTopology, topology(tmesh))
  ntriangles = nelements(t)

  # add middle points of edges
  middles = Dict{Tuple{Int,Int},Int}()
  ∂₁₀ = Boundary{1,0}(t)
  for eind in 1:nfacets(t)
    e1, e2 = ∂₁₀(eind)
    p1, p2 = points[e1], points[e2]
    midpoint = Point((coordinates(p1) + coordinates(p2))/2)
    push!(points, midpoint)
    npoints = npoints + 1
    middles[_ordered_pair(e1, e2)] = npoints
  end

  # construct subtriangles of faces
  ∂₂₀ = Boundary{2,0}(t)
  triangles = Vector{Tuple{Int,Int,Int}}(undef, 4 * ntriangles)
  for tind in 1:ntriangles
    v1, v2, v3 = ∂₂₀(tind)
    m12 = middles[_ordered_pair(v1, v2)]
    m23 = middles[_ordered_pair(v2, v3)]
    m31 = middles[_ordered_pair(v3, v1)]
    k = 4 * (tind - 1)
    triangles[k + 1] = (v1, m12, m31)
    triangles[k + 2] = (v2, m23, m12)
    triangles[k + 3] = (v3, m31, m23)
    triangles[k + 4] = (m12, m23, m31)
  end

  SimpleMesh(points, connect.(triangles, Triangle))
end

_ordered_pair(i ,j) = i < j ? (i, j) : (j, i)

