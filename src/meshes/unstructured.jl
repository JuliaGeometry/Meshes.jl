# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UnstructuredMesh(points, connec)

An unstructured mesh with `points` and connectivities
`connec`. The i-th element of the mesh is lazily built
based on the connectivity list `connec[i]`.
"""
struct UnstructuredMesh{Dim,T,Connectivity} <: Mesh{Dim,T}
  points::Vector{Point{Dim,T}}
  connec::Vector{Connectivity}
end

vertices(m::UnstructuredMesh) = m.points

function faces(m::UnstructuredMesh{Dim,T}, r) where {Dim,T}
  @assert r ≤ Dim "invalid rank for mesh"
  ps, cs = m.points, m.connec
  facerank(c) = rank(facetype(c){Dim,T})
  (materialize(cs[i], ps) for i in eachindex(cs) if facerank(cs[i]) == r)
end
