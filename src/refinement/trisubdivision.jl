# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  TriSubdivision()

Refinement of a mesh by preliminarly triangulating it if needed and
then subdividing each triangle into four triangles.

## References

* Charles Loop. 1987. [Smooth subdivision surfaces based on 
  triangles](https://charlesloop.com/thesis.pdf). 
  Master's thesis, University of Utah.
"""
struct TriSubdivision <: RefinementMethod end

function refine(mesh, ::TriSubdivision)
  assertion(paramdim(mesh) == 2, "TriSubdivision only defined for surface meshes")

  # triangulate mesh if necessary
  tmesh = eltype(mesh) <: Triangle ? mesh : simplexify(mesh)

  # retrieve initial list of points
  points = vertices(tmesh)

  # initial number of points and triangles
  np = nvertices(tmesh)
  nt = nelements(tmesh)

  # convert topology to half-edge structure
  t = convert(HalfEdgeTopology, topology(tmesh))

  # add midpoints of edges
  midpoints = Dict{Tuple{Int,Int},Int}()
  ∂₁₀ = Boundary{1,0}(t)
  for eind in 1:nfacets(t)
    i, j = sort(∂₁₀(eind))
    edge = Segment(points[i], points[j])
    push!(points, center(edge))
    midpoints[(i, j)] = (np += 1)
  end

  # construct subtriangles of faces
  ∂₂₀ = Boundary{2,0}(t)
  triangles = Vector{Tuple{Int,Int,Int}}(undef, 4nt)
  for tind in 1:nt
    t = 4 * (tind - 1)
    i, j, k = ∂₂₀(tind)
    m1 = midpoints[_ordered(i, j)]
    m2 = midpoints[_ordered(j, k)]
    m3 = midpoints[_ordered(k, i)]
    triangles[t + 1] = (i, m1, m3)
    triangles[t + 2] = (j, m2, m1)
    triangles[t + 3] = (k, m3, m2)
    triangles[t + 4] = (m1, m2, m3)
  end

  SimpleMesh(points, connect.(triangles))
end

_ordered(i, j) = i < j ? (i, j) : (j, i)
