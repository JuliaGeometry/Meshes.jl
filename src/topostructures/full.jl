# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FullStructure(connectivities)

A data structure that stores *all* `connectivities` of a mesh.

### Notes

This data structure is sometimes referred to as the "soup of geometries".
It does *not* support topological relations and is therefore incompatible
with algorithms that rely on neighborhood search. It is still useful for
mesh visualization and IO operations.
"""
struct FullStructure{C<:Connectivity} <: TopologicalStructure
  # input fields
  connec::Vector{C}

  # state fields
  ranks::Vector{Int}
  elms::Vector{Int}

  function FullStructure{C}(connec) where {C}
    ranks = [paramdim(c) for c in connec]
    elms  = findall(isequal(maximum(ranks)), ranks)
    new(connec, ranks, elms)
  end
end

FullStructure(connec) = FullStructure{eltype(connec)}(connec)

==(s1::FullStructure, s2::FullStructure) = s1.connec == s2.connec

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(s::FullStructure) = maximum(i for c in s.connec for i in indices(c))

function faces(s::FullStructure, rank)
  cs, rs = s.connec, s.ranks
  (cs[i] for i in 1:length(cs) if paramdim(cs[i]) == rank)
end

element(s::FullStructure, ind) = s.connec[s.elms[ind]]

nelements(s::FullStructure) = length(s.elms)

facets(s::FullStructure) = faces(s, maximum(s.ranks) - 1)

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{FullStructure}, s::TopologicalStructure)
  # TODO: add all faces, not just the elements
  FullStructure(collect(elements(s)))
end
