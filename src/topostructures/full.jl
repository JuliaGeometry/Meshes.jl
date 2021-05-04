# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FullStructure(connectivities)

A data structure that stores *all* `connectivities` of a mesh.
Unlike the [`ElementListStructure`](@ref), this data structure
also stores connectivities of faces that are not elements. For
example, in a mesh of triangle elements this data structure also
stores the connectivities of segments.
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

==(s1::FullStructure, s2::FullStructure) =
  s1.connec == s2.connec

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

function faces(s::FullStructure, rank)
  cs, rs = s.connec, s.ranks
  (cs[i] for i in 1:length(cs) if paramdim(cs[i]) == rank)
end

element(s::FullStructure, ind) = s.connec[s.elms[ind]]
nelements(s::FullStructure) = length(s.elms)
