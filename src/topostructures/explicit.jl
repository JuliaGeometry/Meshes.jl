# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExplicitStructure(connectivities)

A data structure that encodess all incidence relations
of a mesh explicitly as a list of `connectivities`.

See also `[`Connectivity`](@ref).
"""
struct ExplicitStructure{C<:Connectivity} <: TopologicalStructure
  connectivities::Vector{C}
end

connectivities(s::ExplicitStructure) = s.connectivities

==(s1::ExplicitStructure, s2::ExplicitStructure) =
  Set(s1.connectivities) == Set(s2.connectivities)
