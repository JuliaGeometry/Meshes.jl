# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FullTopology(connectivities)

A data structure that stores *all* `connectivities` of a mesh.

### Notes

This data structure is sometimes referred to as the "soup of geometries".
It does *not* support topological relations and is therefore incompatible
with algorithms that rely on neighborhood search. It is still useful for
mesh visualization and IO operations.
"""
struct FullTopology{C<:Connectivity} <: Topology
  # input fields
  connec::Vector{C}

  # state fields
  ranks::Vector{Int}
  elms::Vector{Int}

  function FullTopology{C}(connec) where {C}
    ranks = [paramdim(c) for c in connec]
    elms  = findall(isequal(maximum(ranks)), ranks)
    new(connec, ranks, elms)
  end
end

FullTopology(connec) = FullTopology{eltype(connec)}(connec)

paramdim(t::FullTopology) = paramdim(t.connec[first(t.elms)])

==(t1::FullTopology, t2::FullTopology) = t1.connec == t2.connec

"""
    connec4elem(t, e)

Return linear indices of vertices of `e`-th element of
the full topology `t`.
"""
connec4elem(t::FullTopology, e) = indices(t.connec[t.elms[e]])

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::FullTopology) = maximum(i for c in t.connec for i in indices(c))

function faces(t::FullTopology, rank)
  cs = t.connec
  (cs[i] for i in 1:length(cs) if paramdim(cs[i]) == rank)
end

element(t::FullTopology, ind) = t.connec[t.elms[ind]]

nelements(t::FullTopology) = length(t.elms)

facets(t::FullTopology) = faces(t, maximum(t.ranks) - 1)

nfacets(t::FullTopology) = count(==(maximum(t.ranks) - 1), t.ranks)

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{FullTopology}, t::Topology)
  # TODO: add all faces, not just the elements
  FullTopology(collect(elements(t)))
end
