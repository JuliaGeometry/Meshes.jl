# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ElementListStructure(connectivities)

A data structure that only stores the elements (or top-faces) of a mesh.
Each element is represented with a [`Connectivity`](@ref) in a vector of
`connectivities`.

### Notes

This data structure is sometimes referred to as the "soup of geometries".
It does *not* support topological relations and is therefore incompatible
with algorithms that rely on neighborhood search. It is still useful for
mesh visualization and IO operations.
"""
struct ElementListStructure{C<:Connectivity} <: TopologicalStructure
  connec::Vector{C}
end

Base.getindex(s::ElementListStructure, ind) = getindex(s.connec, ind)
Base.length(s::ElementListStructure) = length(s.connec)
Base.eltype(s::ElementListStructure) = eltype(s.connec)
Base.firstindex(s::ElementListStructure) = firstindex(s.connec)
Base.lastindex(s::ElementListStructure) = lastindex(s.connec)
Base.iterate(s::ElementListStructure, state=firstindex(s)) =
  state > length(s) ? nothing : (s[state], state+1)

==(s1::ElementListStructure, s2::ElementListStructure) =
  s1.connec == s2.connec
