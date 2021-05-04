# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ElementStructure(connectivities)

A data structure that only stores the elements (or top-faces) of a mesh.
Each element is represented with a [`Connectivity`](@ref) in a vector of
`connectivities`.

### Notes

This data structure is sometimes referred to as the "soup of geometries".
It does *not* support topological relations and is therefore incompatible
with algorithms that rely on neighborhood search. It is still useful for
mesh visualization and IO operations.
"""
struct ElementStructure{C<:Connectivity} <: TopologicalStructure
  connec::Vector{C}
end

==(s1::ElementStructure, s2::ElementStructure) =
  s1.connec == s2.connec

# ---------------------
# HIGH-LEVEL INTERFACE
# ---------------------

element(s::ElementStructure, ind) = s.connec[ind]
nelements(s::ElementStructure) = length(s.connec)
