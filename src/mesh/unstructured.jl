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

function faces(m::UnstructuredMesh{Dim}, r) where {Dim}
  @assert 0 ≤ r ≤ Dim "invalid rank for mesh"
  ps, cs = m.points, m.connec
  r == 0 && return ps
  rank(c) = paramdim(polytopetype(c))
  (materialize(c, ps) for c in cs if rank(c) == r)
end

function Base.show(io::IO, m::UnstructuredMesh{Dim,T}) where {Dim,T}
  nvert = length(m.points)
  nface = length(m.connec)
  print(io, "UnstructuredMesh($nvert vertices, $nface faces)")
end

function Base.show(io::IO, ::MIME"text/plain", m::UnstructuredMesh{Dim,T}) where {Dim,T}
  nvert = length(m.points)
  nface = length(m.connec)
  println(io, "UnstructuredMesh{$Dim,$T}")
  println(io, "  $nvert vertices")
  lines = ["    └─$p" for p in m.points]
  println(io, join(lines, "\n"))
  println(io, "  $nface faces")
  lines = ["    └─$f" for f in m.connec]
  print(  io, join(lines, "\n"))
end
