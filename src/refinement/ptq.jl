# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PTQ()

Refinement of triangle meshes.
A triangle is subdivided into four triangles.
"""
struct PTQ <: RefinementMethod end

function refine(mesh, ::PTQ)
  primal_triangle_quadrisection(mesh, (p1, p2) -> Point((coordinates(p1) + coordinates(p2))/2))
end

function ordered_pair(i ,j)
  i < j ? (i, j) : (j, i)
end

function primal_triangle_quadrisection(mesh, middle)
  if !(eltype(mesh) <: Triangle)
    throw(ArgumentError("The mesh is not triangle."))
  end

  # retrieve geometry and topology
  points = vertices(mesh)
  npoints = length(points)
  connec = topology(mesh)

  # convert to half-edge structure
  t = convert(HalfEdgeTopology, connec)
  ntriangles = nelements(t)

  # add middle points of edges
  middles = Dict{Tuple{Int,Int},Int}()
  ∂₁₀ = Boundary{1,0}(t)
  for edge_index in 1:nfacets(t)
    e1, e2 = ∂₁₀(edge_index)
    push!(points, middle(points[e1], points[e2]))
    npoints = npoints + 1
    middles[ordered_pair(e1, e2)] = npoints
  end

  # construct subtriangles of faces
  ∂₂₀ = Boundary{2,0}(t)
  triangles = Vector{Tuple{Int,Int,Int}}(undef, 4 * ntriangles)
  for face_index in 1:ntriangles
    v1, v2, v3 = ∂₂₀(face_index)
    m12 = middles[ordered_pair(v1, v2)]
    m23 = middles[ordered_pair(v2, v3)]
    m31 = middles[ordered_pair(v3, v1)]
    k = 4 * (face_index - 1)
    triangles[k + 1] = (v1, m12, m31)
    triangles[k + 2] = (v2, m23, m12)
    triangles[k + 3] = (v3, m31, m23)
    triangles[k + 4] = (m12, m23, m31)
  end

  SimpleMesh(points, connect.(triangles, Triangle))
end
