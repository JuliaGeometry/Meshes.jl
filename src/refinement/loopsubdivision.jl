# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  LoopSubdivision()

Refinement of a mesh by preliminarly triangulating it if needed and
then subdividing each triangle into four triangles.

## References

* Charles Loop. 1987. [Smooth subdivision surfaces based on 
  triangles](https://charlesloop.com/thesis.pdf). 
  Master's thesis, University of Utah.
"""
struct LoopSubdivision <: RefinementMethod end

function refine(mesh, ::LoopSubdivision)
  tmesh = eltype(mesh) <: Triangle ? mesh : simplexify(mesh)

  # retrieve vertices
  points = vertices(tmesh)
  npoints = length(points)

  # convert topology to half-edge structure
  t = convert(HalfEdgeTopology, topology(tmesh))
  ntriangles = nelements(t)

  # add midpoints of edges
  midpoints = Dict{Tuple{Int,Int},Int}()
  ∂₁₀ = Boundary{1,0}(t)
  for eind in 1:nfacets(t)
    i, j = sort(∂₁₀(eind))
    midpoint = Segment(points[i], points[j])(0.5)
    push!(points, midpoint)
    npoints = npoints + 1
    midpoints[(i, j)] = npoints
  end

  # construct subtriangles of faces
  ∂₂₀ = Boundary{2,0}(t)
  triangles = Vector{Tuple{Int,Int,Int}}(undef, 4 * ntriangles)
  for tind in 1:ntriangles
    v1, v2, v3 = ∂₂₀(tind)
    m12 = midpoints[_ordered_pair(v1, v2)]
    m23 = midpoints[_ordered_pair(v2, v3)]
    m31 = midpoints[_ordered_pair(v3, v1)]
    k = 4 * (tind - 1)
    triangles[k + 1] = (v1, m12, m31)
    triangles[k + 2] = (v2, m23, m12)
    triangles[k + 3] = (v3, m31, m23)
    triangles[k + 4] = (m12, m23, m31)
  end

  SimpleMesh(points, connect.(triangles, Triangle))
end

_ordered_pair(i ,j) = i < j ? (i, j) : (j, i)

