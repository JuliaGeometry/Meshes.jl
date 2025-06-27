# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleTopology(connectivities)

A data structure that stores *all* `connectivities` of a mesh.

### Notes

This data structure is sometimes referred to as the "soup of geometries".
It does *not* support topological relations and is therefore incompatible
with algorithms that rely on neighborhood search. It is still useful for
mesh visualization and IO operations.
"""
struct SimpleTopology{C<:Connectivity} <: Topology
  # input fields
  connec::Vector{C}

  # state fields
  ranks::Vector{Int}
  elems::Vector{Int}

  function SimpleTopology{C}(connec) where {C}
    if isconcretetype(C) # single concrete connectivity type
      ranks = fill(paramdim(first(connec)), length(connec))
      elems = collect(eachindex(connec))
    else # mixed connectivity types
      ranks = [paramdim(c) for c in connec]
      elems = findall(isequal(maximum(ranks)), ranks)
    end
    new(connec, ranks, elems)
  end
end

SimpleTopology(connec) = SimpleTopology{eltype(connec)}(connec)

paramdim(t::SimpleTopology) = paramdim(t.connec[first(t.elems)])

==(t1::SimpleTopology, t2::SimpleTopology) = t1.connec == t2.connec

"""
    connec4elem(t, e)

Return linear indices of vertices of `e`-th element of
the simple topology `t`.
"""
connec4elem(t::SimpleTopology, e) = indices(t.connec[t.elems[e]])

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

nvertices(t::SimpleTopology) = maximum(i for c in t.connec for i in indices(c))

function faces(t::SimpleTopology, rank)
  cs = t.connec
  (cs[i] for i in 1:length(cs) if paramdim(cs[i]) == rank)
end

element(t::SimpleTopology, ind) = t.connec[t.elems[ind]]

nelements(t::SimpleTopology) = length(t.elems)

facets(t::SimpleTopology) = faces(t, maximum(t.ranks) - 1)

nfacets(t::SimpleTopology) = count(==(maximum(t.ranks) - 1), t.ranks)

# ------------
# CONVERSIONS
# ------------

function Base.convert(::Type{SimpleTopology}, t::Topology)
  ranksₜ = 1:paramdim(t)
  facesₜ(r) = collect(faces(t, r))
  connec = mapreduce(facesₜ, vcat, reverse(ranksₜ))
  SimpleTopology(connec)
end
