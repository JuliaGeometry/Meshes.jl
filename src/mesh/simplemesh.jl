# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, connec)

A simple mesh with `points` and connectivities `connec`.
The i-th face of the mesh is lazily built based on
the connectivity list `connec[i]`.
"""
struct SimpleMesh{Dim,T,Connectivity} <: Mesh{Dim,T}
  # input fields
  points::Vector{Point{Dim,T}}
  connec::Vector{Connectivity}

  # state fields
  ranks::Vector{Int}
  elms::Vector{Int}

  function SimpleMesh{Dim,T,C}(points, connec) where {Dim,T,C}
    ranks = [paramdim(polytopetype(c)) for c in connec]
    elms  = [i for i in eachindex(ranks) if ranks[i] == Dim]
    new(points, connec, ranks, elms)
  end
end

function SimpleMesh(points, connec)
  p = first(points)
  Dim = embeddim(p)
  T = coordtype(p)
  C = eltype(connec)
  SimpleMesh{Dim,T,C}(points, connec)
end

==(m1::SimpleMesh, m2::SimpleMesh) =
  m1.points == m2.points && m1.connec == m2.connec

vertices(m::SimpleMesh) = m.points

function faces(m::SimpleMesh{Dim}, r) where {Dim}
  @assert 0 < r ≤ Dim "invalid rank for mesh"
  ps, cs, rs = m.points, m.connec, m.ranks
  (materialize(cs[i], ps) for i in eachindex(cs) if rs[i] == r)
end

# -----------------
# DOMAIN INTERFACE
# -----------------

Base.getindex(m::SimpleMesh, ind::Int) =
  materialize(m.connec[m.elms[ind]], m.points)

nelements(m::SimpleMesh) = length(m.elms)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::MIME"text/plain", m::SimpleMesh{Dim,T}) where {Dim,T}
  nvert = length(m.points)
  nface = length(m.connec)
  println(io, m)
  println(io, "  $nvert vertices")
  println(io, _lines(m.points, "    "))
  println(io, "  $nface faces")
  print(  io, _lines(m.connec, "    "))
end

function _lines(vec, tab="  ")
  N = length(vec)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["$(tab)└─$(vec[i])" for i in 1:I]
           (N > 10 ? ["$(tab)⋮"] : [])
           ["$(tab)└─$(vec[i])" for i in J:N]]
  join(lines, "\n")
end
