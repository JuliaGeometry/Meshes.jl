# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
Convention used to reconstruct the faces of a polytope
from a list of vertices.
"""
abstract type OrderingConvention end

"""
    connectivity(polytope, convention, k)

Return the list of `Connectivity` elements used to construct all
faces of rank `rank` of a polytope, following an ordering convention.
"""
function connectivity(polytope::Type{<:Polytope}, rank::Integer, convention::OrderingConvention)
  connectivity(polytope, Val(rank), convention)
end

# ---------------------------
# CONVENTION IMPLEMENTATIONS
# ---------------------------

include("conventions/gmsh.jl")

function facets(p::Polytope{N}, ordering=GMSH) where {N}
  faces(p, N-1, ordering)
end

function faces(p::Polytope, rank, ordering=GMSH)
  (materialize(ord, p.vertices) for ord in connectivity(typeof(p), rank, ordering))
end
