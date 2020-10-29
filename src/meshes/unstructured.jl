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

function Base.show(io::IO, m::UnstructuredMesh{Dim,T}) where {Dim,T}
  nvert = length(m.points)
  nface = length(m.connec)
  print(io, "UnstructuredMesh($nvert vertices, $nface faces)")
end

function Base.show(io::IO, ::MIME"text/plain", m::UnstructuredMesh{Dim,T}) where {Dim,T}
  println(io, "UnstructuredMesh{$Dim,$T}")
  for r in 1:Dim
    fs = collect(faces(m, r))
    n  = length(fs)
    if n > 0
      lines = ["    └─$f" for f in fs]
      println(io, "  $r-faces")
      print(  io, join(lines, "\n"))
    end
  end
end
